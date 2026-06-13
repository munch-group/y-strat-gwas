# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What this is

A two-arm genetic analysis pipeline (a `gwf` workflow + small Python helpers) testing for
**SNP × Y-haplogroup epistasis** in a male autism cohort split into Y-haplogroups **I** and **R**.
The deliverable is `workflow.py` plus the helper scripts it invokes, executed on the GenomeDK
Slurm cluster. `tests/` holds a self-contained local smoke test (synthetic data, no Slurm).

`METHODS.md` is the standalone conceptual guide (the statistics and *why* the pipeline is built
this way) — read it for the theory; this file is the operational map. The single most important
domain fact (see `METHODS.md`, `notes.md`, README caveat 1): haplogroup I vs R is
**not exogenous** — it tags deep autosomal ancestry (R carries more Steppe-related ancestry).
Every design choice exists to separate genuine epistasis from this ancestry confounding. When
editing, preserve those guards rather than "simplifying" them away:
- The `PC_i × Hap` product covariates in Arm A (Keller 2014 E×covariate guard).
- The warning against case-only designs.
- `STRATUM_SPECIFIC_STEP1` and its r_g-toward-1 bias note.

## Two arms

- **Arm A — individual variants** (`step1_full` → `interaction` → `top_int`): full-sample
  REGENIE step-1, then a genome-wide step-2 `--interaction Hap` scan whose covariate file
  carries `Hap` + `PC1_x_Hap … PCk_x_Hap`. `top_interactions.py` extracts the `ADD-INT_SNPxHap`
  rows and reports λ_GC.
- **Arm B — general pattern** (`gwas_{I,R}` → `munge_{I,R}` → `h2_{I,R}`, plus `rg`): a
  main-effect GWAS within each stratum, munged to LDSC format, then cross-stratum genetic
  correlation. **r_g significantly < 1 ⇒ pervasive G × haplogroup interaction** — the
  well-powered readout. Arm B optionally runs a stratum-specific step-1 (`STRATUM_SPECIFIC_STEP1`).

Coding convention enforced across all helpers: **Hap is I=1, R=0**. `ALLELE1`/`ALLELE0` are
REGENIE's effect/non-effect alleles; REGENIE reports `LOG10P`, not `P` (helpers convert).

- **Ancestry-matched permutation** (`perm_interaction`, post-hoc to Arm A): the G-side guard
  REGENIE's `--interaction` can't add. It re-tests interaction hits against a null that permutes
  Hap **within Hap-propensity strata** — a real I-vs-R effect lives within strata and is broken
  (survives; small `anc_matched_emp_p`), an ancestry artifact lives between strata and is
  preserved (collapses; p→1). Strata granularity is auto-selected: fine enough that Hap is
  unpredictable from PCs *within* a stratum (AUC→0.5), then the loosest such for power. Crucially
  it controls confounding only to **global-PC resolution** — local ancestry at a locus survives
  it, so survivors are still provisional. Stratifying on the 1-D propensity (not raw multi-PC
  distance) is what makes decoupling reachable.
- **LAVA** (optional, only when `LAVA_LOCI` is set): local cross-stratum r_g per LD block —
  the local-resolution version of Arm B, for following up candidate regions. Same ancestry
  confound (differential LD depresses local r_g genome-wide), so the headline LAVA estimate
  (`lava_local`, the real **LAVA R package** — not in the pixi env, handled like LDSC) is paired
  with a pure-Python deconfounding layer (`local_rg_perm.py`): the Section-6 permutation applied
  to an LD-aware cross-stratum beta-correlation, plus **region-excluded PCs** (anti-circularity)
  and a **negative-control-loci panel** (a focused locus is only interesting if it's in the tail
  of matched control blocks, since residual structure depresses local r_g everywhere). Same
  global-PC ceiling — local ancestry survives.
- **chrX** (optional, only when `XBFILE` is set): both arms get a chrX step-2 pass reusing the
  autosomal step-1 predictors, with non-PAR males coded **hemizygously (0/2)** by REGENIE from
  the `.fam` sex column + `--par-region <build>` (`--sex-specific male`). Targets:
  `xqc` → `interaction_X` → `top_int_X`, and per stratum `gwas_X_{I,R}` → `xforldsc_{I,R}`.
  LDSC stays autosomal — `rg` does not include chrX; the `gwas_X_{I,R}.forldsc.txt` sumstats are
  emitted munge-ready, but X h2/rg needs an X-specific LD-score reference the user supplies.

## Single pixi environment

Everything — regenie, plink2, pandas, gwf, the helper scripts, **and LDSC** — runs under one
pixi env (`pixi.toml`, Python 3), invoked as `pixi run --manifest-path .../pixi.toml`. LDSC is a
**Python 3 fork vendored in `./ldsc`** and pip-installed into the pixi env by the `postinstall`
task (`pip install --no-deps --no-build-isolation ./ldsc`); its C-extension deps (bitarray,
pybedtools) come from conda so pip never builds anything. `LDSCRUN` is therefore just
`<pixi run> python` and `LDSC_DIR` defaults to `./ldsc`. The legacy `environment-ldsc.yml`
(Python 2.7) is **not used** by the workflow and is kept only as a fallback reference.

The one exception is the optional **LAVA R headline** (`lava_local` / `lava_local.R`): LAVA is
an R package, not conda-native, so it runs via `RSCRIPT` outside the pixi env (install it
separately, like a vendored dependency). Everything else in the LAVA arm — the deconfounding,
permutation and negative controls (`local_rg_perm.py`) — is pure Python in the pixi env.

## Setup

```bash
pixi install   # all tools; postinstall pip-installs the vendored ./ldsc fork into the env
# LDSC repo is vendored in ./ldsc (committed without git history); LDSC_DIR defaults to it.
```

LDSC reference data must be downloaded separately and pointed at via `workflow.py`:
`eur_w_ld_chr/` → `EUR_LD`, `w_hm3.snplist` → `HM3`.

## Configure & run

All paths and parameters live in the **CONFIG block at the top of `workflow.py`**
(`ACCOUNT`, `BFILE`, `PHENO`, `BASECOVAR`, `HAPFILE`, `NPC`, `CATCOVAR`, `LDSC_DIR`, `EUR_LD`,
`HM3`, `PREV_POP`, `STRATUM_SPECIFIC_STEP1`, `XBFILE`/`GENOME_BUILD` for chrX, `PERM_*` for the
permutation re-test, and `LAVA_*`/`RSCRIPT` for the LAVA arm). These ship as `/path/to/...`
placeholders (`XBFILE` and `LAVA_LOCI` default empty = those arms skipped). Each can
be **edited in place or overridden via a `YS_<NAME>` environment variable** (e.g. `YS_BFILE`,
`YS_NPC`, `YS_OUT`) — the env path is what the test harness uses, and it avoids editing the file
on the cluster.

```bash
gwf -f workflow.py status
gwf -f workflow.py run
```

Outputs land in `results/`; REGENIE step-1 lowmem scratch in `tmp/`. Key results:
`results/top_interactions.tsv`, `results/I_vs_R_rg.log`, `results/h2_{I,R}.log`.

The helper scripts live in `scripts/` and `workflow.py` invokes them as `{ROOT}/scripts/<name>.py`.

## Helper scripts (each is a standalone CLI, run via pixi)

| script | role |
|--------|------|
| `make_strata.py` | Hap file → per-stratum `keep_{I,R}.txt` for REGENIE `--keep` |
| `make_interaction_covars.py` | base covars + Hap → `int_covars.txt` with `PC×Hap` products |
| `regenie_to_munge.py` | REGENIE step-2 → LDSC-munge-ready table (LOG10P→P, ALLELE1→A1) |
| `top_interactions.py` | extract `ADD-INT_SNPxHap` rows, report λ_GC, dump top hits |
| `parse_ldsc_rg.py` | print the genetic-correlation summary block from an `--rg` log |
| `samp_prev.py` | per-stratum case fraction → LDSC `--samp-prev` (required with `--pop-prev`) |
| `ancestry_matched_perm.py` | re-test interaction hits vs Hap-within-propensity-strata permutation null; reads `.bed` directly |
| `local_rg_perm.py` | LAVA local-r_g deconfounding: ancestry-matched permutation + negative-control loci (imports `ancestry_matched_perm`) |
| `lava_inputs.py` | build LAVA input.info / sample-overlap / locus files from the per-stratum sumstats |
| `lava_local.R` | headline bivariate local r_g via the LAVA R package (run outside the pixi env) |

`top_interactions.py --test` defaults to `ADD-INT_SNPxHap`; **verify the exact TEST label in
your REGENIE output header** — it varies by version, and a mismatch makes the script exit empty.

## Local test (`tests/`)

`pixi run --manifest-path ./pixi.toml python tests/run_pipeline_test.py` runs the **whole
workflow** on tiny synthetic data with no Slurm: it generates a confounded I/R dataset with one
planted SNP×Hap interaction, points `workflow.py` at it via `YS_*` env vars, executes each gwf
target's `spec` in order, and asserts the planted SNP lands in the top 5% of the scan plus that
rg / liability-h2 come out. The synthetic data also plants a **confound-only** SNP (effect varies
with ancestry, no true Hap interaction) and the test asserts the ancestry-matched permutation
**separates** it from the real interaction (real survives, confound collapses). The data also
plants two LD blocks — a **genuine** local divergence (effects sign-flip by Hap) and a **confound**
one (effects modulated by ancestry) — and the test asserts the LAVA permutation + negative-control
layer separates them (genuine survives and beats controls; confound collapses). It also builds a
synthetic chrX fileset (all-male, hemizygous) and checks the planted **chrX** interaction surfaces
from the hemizygous step-2 pass. The `lava_local` R step is **not** exercised (no R here). Scratch goes to
`tests/work/` (git-ignored). See `tests/README.md`. Use this to smoke-test any change to
`workflow.py` or the helpers before a cluster run.

The vendored `./ldsc` fork needed small patches to run under modern pandas/numpy (all genuine
py2→py3 / pandas-2 / numpy-2 breakages that bite real runs, found via the test):
- `ldscore/sumstats.py`: `.drop('SNP', 1)` → `.drop('SNP', axis=1)` (pandas 2 positional axis).
- `ldscore/regressions.py`: `float(rg.jknife_est)` etc. → `float(np.asarray(x).ravel()[0])`
  (numpy 2 won't convert a size-1 `(1,1)` array to a scalar — the `rg` path hits this whenever
  the h2 product is positive, i.e. the normal case).
- `ldscore/sumstats.py`: `traceback.format_exc(ex)` → `traceback.format_exc()` (py3 signature;
  the old call crashed the rg error-handler with a misleading secondary `TypeError`).

Its bed reader (`--l2`) is still incompatible with the installed `bitarray`, but the workflow
never calls `--l2` (it consumes pre-made reference LD scores), so that path is unused.

## Known stubs / deliberately out of scope

- **chrX LDSC.** The chrX REGENIE step-2 passes exist (hemizygous coding), but LDSC h2/rg on
  chrX is not wired — it needs an X-chromosome LD-score reference + X-inclusive merge-alleles
  list. The `gwas_X_{I,R}.forldsc.txt` outputs are emitted ready to feed once that's available.
- The binary I/R label is a coarse proxy for the real hypothesis (Y heterochromatin dosage).
  A continuous Yq12 read-depth measure could be swapped in wherever `Hap` appears.
- A two-step interaction screen (marginal-variance/association filter → interaction on survivors)
  is discussed but not implemented; the genome-wide pass is underpowered without it.
