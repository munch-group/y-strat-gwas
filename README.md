# GxHaplogroup: SNP × Y-haplogroup interaction + cross-stratum genetic correlation

Two analyses from one plink fileset of male autism cases/controls split into
Y-haplogroup **I** and **R**:

- **Arm A — individual variants.** Genome-wide SNP × haplogroup interaction scan
  in REGENIE, with `PC_i × Hap` product terms included so that autosomal
  ancestry (which differs between I and R lineages — e.g. more Steppe-related
  ancestry in R) cannot masquerade as interaction.
- **Arm B — general pattern.** Main-effect GWAS within each stratum, then LDSC
  cross-stratum genetic correlation. **r_g significantly below 1 ⇒ pervasive
  G × haplogroup interaction.** Per-stratum SNP-h2 is reported too.
- **LAVA (optional) — local pattern.** Local cross-stratum genetic correlation per
  LD block (the local-resolution version of Arm B), paired with an ancestry-matched
  permutation test + negative-control loci to deconfound focused per-locus claims.

> **New here?** Read [`METHODS.md`](METHODS.md) first — a standalone guide to the
> statistical ideas (ancestry confounding, Keller's interaction terms, cross-stratum
> r_g, the ancestry-matched permutation guard) and why the pipeline is built this way.
> This README covers how to *run* it; `METHODS.md` covers *why*.

## Input formats

| file | columns |
|------|---------|
| `genotypes.{bed,bim,fam}` | plink1, autosomes, QC'd samples |
| `phenotypes.txt` | `FID IID autism` — `autism` = 1 case / 0 control / `NA` |
| `covariates_base.txt` | `FID IID PC1 … PCk age batch …` |
| `haplogroup.txt` | `FID IID Hap` — `Hap` = `I` or `R` |
| `genotypesX.{bed,bim,fam}` | *(optional)* plink1 chrX (chromosome coded 23/X), **SEX filled in the `.fam`** (1=male) |

`Hap` is coded internally as I=1, R=0. The optional chrX fileset (`XBFILE`) is analysed in
its own REGENIE step-2 pass — REGENIE codes non-PAR males hemizygously (0/2) from the `.fam`
sex column plus the `--par-region` build, so the sex column **must** be set.

## Setup

```bash
# all tools, including the python3 LDSC fork
pixi install                          # reads pixi.toml (regenie, plink2, py3, gwf)
                                      # postinstall pip-installs the vendored ./ldsc fork
```

LDSC is a python3 fork vendored in `./ldsc` and pip-installed into the pixi env by the
`postinstall` task, so it runs under the same `pixi run` as everything else — no separate
python 2.7 environment is needed. `LDSC_DIR` defaults to `./ldsc` in `workflow.py`.

Download LDSC reference data (distributed by the LDSC authors; mirrored on
Zenodo) and set the paths in `workflow.py`:

- `eur_w_ld_chr/` — 1000G EUR LD scores → `EUR_LD`
- `w_hm3.snplist` — HapMap3 SNPs for `--merge-alleles` → `HM3`

The **LAVA arm is optional** (disabled unless `LAVA_LOCI` is set). Its permutation /
negative-control deconfounding (`lava_perm_*`) is pure Python and runs in the pixi env.
The *headline* local r_g (`lava_local`) uses the **LAVA R package**, which is not in the
pixi env — install it separately (e.g. `R -e 'remotes::install_github("josefin-werme/LAVA")'`),
point `RSCRIPT` at that R, and supply a `LAVA_PARTITION` blocks file; only then is it built.

## Configure & run

The CONFIG block at the top of `workflow.py` **defaults to the bundled synthetic test
data** (`tests/work/data`) so the workflow runs out of the box once that data is built
(see the smoke test below); `NPC=10` and `ENV_PREFIX=KMP_AFFINITY=disabled` are set to
match. **Edit the CONFIG block to point at your real data** (`BFILE`, `PHENO`,
`BASECOVAR`, `HAPFILE`, `EUR_LD`, `HM3`, `NPC`, `PREV_POP`, `ACCOUNT`, …) and clear
`ENV_PREFIX` for a real cluster run. To enable the chrX analysis set `XBFILE` (and
`GENOME_BUILD`, default `hg38`); to enable the LAVA arm set `LAVA_LOCI`
(e.g. `HCN1=5:45300000-45700000`). Both are empty (skipped) by default.

Any CONFIG value can instead be supplied via a `YS_<NAME>` environment variable (e.g.
`YS_BFILE`, `YS_NPC`), so you don't have to edit the file on the cluster. `ENV_PREFIX`
prepends extra environment to every command — e.g.
`YS_ENV_PREFIX="KMP_AFFINITY=disabled"` to dodge an MKL/OpenMP thread-affinity assertion
in some sandboxes/containers; because it's baked into the spec strings it also applies
under a `gwf` worker whose own environment lacks it.

```bash
gwf -f workflow.py status
gwf -f workflow.py run
```

### Parallelising on a cluster (optional)

The strata, arms and LAVA loci already run as independent `gwf` jobs. Two config flags add
data-subset parallelism (both **on by default**; both leave results unchanged). Set
`SPLIT_CHROMS=""` and `PERM_BATCHES=1` for the simple one-job-per-scan path.

- **`SPLIT_CHROMS`** (default `"1-22"`; also accepts e.g. `"1,2,5-9"`) — fans the REGENIE step-2
  scans (`interaction`, `gwas_{I,R}`) into one job per chromosome plus a gather. Step 2 tests
  variants independently, so the gathered output is **identical** to the single job (verified in
  `tests/test_parallel.py`).
- **`PERM_BATCHES`** (default `4`) — splits the Monte-Carlo permutation targets (`perm_interaction`,
  `lava_perm_*`) into N batches with distinct permutation seeds (selection kept fixed), pooled
  into the same final result. Linear speedup for large `PERM_NPERM`/`LAVA_NPERM`.

On the synthetic data (chr 1–22, both optional arms on) this fans out to ~104 `gwf` tasks; the
whole split pipeline is verified to complete with correct gathered/pooled results.

(REGENIE step 1 is a joint fit and doesn't slice by chromosome; use REGENIE's own
`--split-l0`/`--run-l0`/`--run-l1` if step 1 ever becomes the bottleneck.)

### Smoke-test the whole pipeline locally (no cluster)

```bash
pixi run --manifest-path ./pixi.toml python tests/run_pipeline_test.py
```

This runs every target's command in-process. To instead exercise the **real `gwf` local
backend** (worker daemon + dependency scheduling — a true cluster emulation) against the
same synthetic data, run `bash tests/run_via_gwf.sh` after the smoke test has built the
data once. Both have been verified to complete all targets and produce correct results.

Runs every stage on tiny synthetic data and checks a planted SNP×Hap interaction surfaces in
the scan and that rg / liability-h2 are produced. See `tests/README.md`.

Target graph: `qc`, `strata`, `int_covars` → `step1_full` →
`interaction` → `top_int` → `perm_interaction`; and per stratum
`gwas_{I,R}` → `munge_{I,R}` → `h2_{I,R}`, `rg`. With `XBFILE` set, the chrX arm
adds `xqc` → `interaction_X` → `top_int_X` (reusing `step1_full`) and per stratum
`gwas_X_{I,R}` → `xforldsc_{I,R}` (reusing the autosomal step-1 predictors). With
`LAVA_LOCI` set, the LAVA arm adds `lava_pcs` → `lava_perm_<name>` (deconfounding) and,
if `LAVA_PARTITION` is also set, `lava_inputs` → `lava_local` (headline local r_g in R).

## Reading the results

See **[`OUTPUTS.md`](OUTPUTS.md)** for a full field-by-field guide to every result file (it
applies to both `results/` and the test's `tests/work/results/`). The key files in brief:

- `results/top_interactions.tsv` — strongest SNP×Hap signals, plus `lambda_GC`
  on the interaction p-values printed to the log. If λ > ~1.10, revisit QC and
  the PC×Hap covars before trusting hits. Genome-wide interaction threshold is
  5e-8 unless you used the two-step screen (below).
- `results/perm_interactions.tsv` — interaction hits re-tested against an
  **ancestry-matched permutation** null (the G-side guard REGENIE can't add).
  `anc_matched_emp_p` is the deconfounded p-value: a real I-vs-R hit stays small,
  an ancestry artifact inflates toward 1. `results/perm_lambda.txt` reports the
  genome-wide interaction λ vs the permutation null, and
  `results/perm_strata_selection.tsv` the strata-granularity choice (see caveat 1).
- `results/I_vs_R_rg.log` — the cross-stratum r_g (parsed to stdout by `rg`).
- `results/h2_{I,R}.log` — per-stratum SNP-h2.
- `results/lava_<name>_summary.txt` + `_loci.tsv` — per-locus **local** cross-stratum
  r_g with its **ancestry-matched** empirical p and its position in the negative-control
  distribution (`target_vs_controls_tail_frac` — small ⇒ stands out from controls). The
  headline LAVA estimate, if you ran the R step, is in `results/lava_local_rg.tsv`.
- `results/top_interactions_X.tsv` — strongest chrX SNP×Hap signals (if `XBFILE` set).
- `results/gwas_X_{I,R}.forldsc.txt` — per-stratum chrX summary stats in LDSC-munge-ready
  form. Running LDSC h2/rg on chrX additionally needs an X-chromosome LD-score reference and
  an X-inclusive `--merge-alleles` list (the bundled autosomal `eur_w_ld_chr` / HapMap3 list
  do not cover chrX), so that final step is left for you to wire up once you have an X reference.

## Important caveats baked into / around this pipeline

1. **Structure is the main threat.** PC×Hap terms guard the E-side; REGENIE
   does not add per-SNP SNP×PC (G-side) terms. The `perm_interaction` target
   supplies that guard post-hoc: it re-tests hits against an **ancestry-matched
   permutation** null (permute Hap within Hap-propensity strata), where a real
   I-vs-R hit survives and an ancestry artifact collapses. **Important:** this
   controls confounding only to **global-PC resolution** — residual *local*
   ancestry at a locus is not removed, so treat survivors as provisional pending
   a local-ancestry re-test or replication. Strata are made fine enough that Hap
   is unpredictable from PCs within a stratum (the AUC criterion in
   `perm_strata_selection.tsv`), then taken as loose as possible for power.
2. **Power.** Interaction needs ~4× the N of a main effect. Consider a two-step
   screen (filter on marginal variance/association, then test interaction only
   in survivors) to cut the multiple-testing burden — or restrict the scan to a
   hypothesis-driven set (imprinted clusters, X-linked regulators) where your
   heterochromatin-sink prior predicts signal. Don't run case-only here: it
   assumes G ⟂ Hap in the population, which the ancestry correlation violates.
3. **r_g standard errors.** Each stratum is ~half the cohort, so LDSC SEs will
   be wide unless N per group is large (rule of thumb: ≳5k per stratum for a
   usable h2/rg). The shared-step1 default biases r_g slightly toward 1
   (conservative); flip `STRATUM_SPECIFIC_STEP1` for the cleaner estimate.
4. **Scale dependence.** Statistical interaction depends on scale (logit vs
   liability). Weight qualitative (sign-flipping) interactions over removable
   quantitative ones.
5. **X chromosome.** Analysed in its own REGENIE step-2 pass (`XBFILE`), with
   non-PAR males coded hemizygously (0/2) via the `.fam` sex column + `--par-region`.
   This covers the interaction scan (`interaction_X`/`top_int_X`) and per-stratum
   chrX association (`gwas_X_{I,R}`). LDSC stays autosomal, so the cross-stratum r_g
   does **not** include chrX — the chrX sumstats are emitted ready to munge, but X
   h2/rg needs an X-specific LD-score reference you supply.
6. **Mechanism.** The binary I/R label is a coarse proxy. A continuous
   Y-heterochromatin estimate (relative read depth over DYZ1/DYZ2) may capture
   the sink effect better — swap it in as the interacting covariate by replacing
   `Hap` with that column.
