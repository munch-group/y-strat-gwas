#!/usr/bin/env bash
# One-command local "cluster" emulation: run the REAL workflow.py end-to-end
# against the synthetic test data via the gwf local backend (a worker daemon),
# exactly as it would run on Slurm but on this machine.
#
#   bash tests/run_via_gwf.sh
#
# It points workflow.py at tests/work/data (via tests/env.sh), starts a gwf worker
# pool, runs all targets, waits for completion, then tears the pool down. The test
# data must already exist -- run tests/run_pipeline_test.py once first (it builds
# the plink filesets, chrX, w_hm3.snplist and the LDSC reference).
#
# For an *interactive* session, `source tests/env.sh` and run gwf yourself.
set -euo pipefail
cd "$(dirname "$0")/.."
ROOT=$PWD
PIXI="pixi run --manifest-path $ROOT/pixi.toml"
PORT=${GWF_PORT:-12345}
NWORK=${GWF_WORKERS:-4}

if [ ! -f "$ROOT/tests/work/data/genotypes.bed" ]; then
  echo "test data not found -- run 'python tests/run_pipeline_test.py' first" >&2
  exit 1
fi
source "$ROOT/tests/env.sh"            # exports all YS_* (incl. YS_ENV_PREFIX=KMP...)
rm -rf "$YS_OUT" "$YS_TMP"; mkdir -p "$YS_OUT" "$YS_TMP"

echo ">> starting $NWORK gwf workers on port $PORT"
$PIXI gwf workers -n "$NWORK" -p "$PORT" >/tmp/gwf_workers.$$.log 2>&1 &
WPID=$!
trap 'kill $WPID 2>/dev/null || true' EXIT
sleep 4

echo ">> submitting all targets"
GWF="$PIXI gwf -b local"
$GWF config set local.port "$PORT" >/dev/null 2>&1 || true
$GWF run >/dev/null 2>&1

echo ">> waiting for completion"
for i in $(seq 1 240); do
  sleep 5
  S=$($GWF status 2>&1 | grep -vE "Could not connect" || true)
  r=$(echo "$S" | grep -c running || true); s=$(echo "$S" | grep -c shouldrun || true)
  f=$(echo "$S" | grep -c failed || true); c=$(echo "$S" | grep -c completed || true)
  printf "  [%3ds] completed=%s running=%s shouldrun=%s failed=%s\n" $((i*5)) "$c" "$r" "$s" "$f"
  if [ "$f" -gt 0 ]; then echo "FAILED targets present -- see .gwf/logs/"; exit 1; fi
  [ "$r" -eq 0 ] && [ "$s" -eq 0 ] && break
  # local backend doesn't always chain a gather after a fast wave finishes;
  # if nothing is running but work remains, re-submit to advance the DAG.
  [ "$r" -eq 0 ] && [ "$s" -gt 0 ] && $GWF run >/dev/null 2>&1 || true
done
echo ">> done. outputs in $YS_OUT"
