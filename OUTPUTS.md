# Reading the output files

This explains the result files the pipeline writes. They land in **`results/`** for a real run
and in **`tests/work/results/`** for the bundled synthetic test (`tests/run_pipeline_test.py`) — the
formats are identical; only the directory differs. For *why* the pipeline is built this way and
how to interpret the statistics, see `METHODS.md`; for what produces each file, the target table
in `CLAUDE.md`.

> **On the synthetic test data the estimates are noise.** N is tiny, so LDSC `h2`/`rg` come out
> negative or `NA` and interaction p-values are not meaningful. The files are still well-formed —
> use them to learn the *format*, not to read off a result. On real data these become informative.

Conventions across every file: **`Hap` is coded I = 1, R = 0**; REGENIE's effect allele is
**`ALLELE1`** (`A1` after munging); REGENIE reports **`LOG10P` = −log₁₀(p)**, not `p`
(`P = 10**(-LOG10P)`; helper scripts already convert where needed).

---

## Start here — `conclusions.txt`

The `summarize` target (runs last) reads every other result file and writes
**`results/conclusions.txt`**: a plain-language statement, per arm, of what the data do and do
**not** support, applying the interpretation rules below (genome-wide significance + permutation
survival for Arm A; the h²-Z / rg-power gate for Arm B; the female-artefact check; LAVA screen +
deconfounding). Every positive is hedged as provisional to the global-PC ceiling. It is a reading
aid, not a substitute for the files it summarises — when a line surprises you, drop to the
underlying file documented below. Disabled or unfinished arms are reported as "not available"
rather than omitted.

---

## Arm A — individual SNP × haplogroup interaction

**`top_interactions.tsv`** — the headline list: the strongest `SNP × Hap` interaction hits.
One row per variant (REGENIE columns), pre-filtered to the interaction test and sorted by
significance:

| column | meaning |
|--------|---------|
| `CHROM`, `GENPOS`, `ID` | variant location and rsID |
| `ALLELE0` / `ALLELE1` | non-effect / **effect** allele |
| `A1FREQ`, `N` | effect-allele frequency, sample size |
| `TEST` | always `ADD-INT_SNPxHap` here (the interaction term) |
| `BETA`, `SE` | interaction effect (on the logit scale) and its SE |
| `CHISQ`, `LOG10P` | test statistic and −log₁₀(p) |
| `P` | the p-value (added by the helper) |

Read it as: large `LOG10P` (small `P`) = strong evidence the variant's effect on autism *differs*
between I and R men. Genome-wide threshold is `5e-8` (i.e. `LOG10P > 7.3`) unless you used a
two-step screen. The genomic-inflation factor **λ_GC** and a warning if λ > 1.10 are printed to
the job's log (the `top_int` stdout), not into the table — λ ≫ 1 means revisit QC / the PC×Hap
covariates before trusting hits.

**`gxhap_autism.regenie`** — the *raw* REGENIE `--interaction` output `top_interactions.tsv` is
distilled from. It has **several `TEST` rows per variant** (`ADD`, `ADD-INT_SNP`,
`ADD-INT_SNPxHap`, …); the interaction you care about is **`ADD-INT_SNPxHap`**. (With chromosome
split on, you may also see per-chunk `gxhap_chr<c>_autism.regenie` — these are intermediate and
already concatenated into `gxhap_autism.regenie`.)

**`perm_interactions.tsv`** — the **deconfounded re-test** of the hits (ancestry-matched
permutation; the G-side guard REGENIE can't add). This is the column that decides whether a hit is
credible:

| column | meaning |
|--------|---------|
| `ID` | variant |
| `raw_LOG10P` | its `LOG10P` from the scan (for reference) |
| `obs_stat` | the observed interaction statistic |
| `n_perm` | permutations used (summed across batches) |
| **`anc_matched_emp_p`** | **empirical p vs the Hap-within-ancestry-strata null** |
| `forced` | whether it was force-included regardless of rank |

**Small `anc_matched_emp_p` ⇒ the interaction survives the ancestry-matched null** (not explained
by *global* ancestry structure). **`anc_matched_emp_p → 1` ⇒ it collapses** = an ancestry
artefact. Caveat (see `METHODS.md` §6): this controls confounding only to *global-PC* resolution;
a survivor is still provisional pending local ancestry or replication.

**`perm_lambda.txt`** — genome-wide companion diagnostic (`key<TAB>value`): the observed
interaction `observed_lambda` vs the permutation-null distribution
(`perm_null_lambda_mean`/`_2.5pct`/`_97.5pct`) and `lambda_emp_p`. A more honest calibration check
than comparing λ to 1. Also records the chosen strata (`selected_stratum_size`,
`within_stratum_HapPC_AUC`, `eff_N`).

**`perm_strata_selection.tsv`** — the strata-granularity choice: for each candidate `target_size`,
the within-stratum `Hap~PC` AUC (want ≈ 0.5), number of mixed strata, and effective N. The run
picks the loosest size whose AUC ≤ threshold (the "tightening" rule).

### chrX (when `XBFILE` is set)
**`top_interactions_X.tsv`** / **`gxhapX_autism.regenie`** — exactly as above but for the
hemizygous-coded chrX pass.

---

## Arm B — genome-wide cross-stratum architecture (LDSC)

**`I_vs_R_rg.log`** — the **cross-stratum genetic correlation**, the well-powered global readout.
Find the `Summary of Genetic Correlation Results` block; the key fields on the data line are
**`rg`**, its **`se`**, and **`p`** (plus per-stratum `h2_obs` and the intercepts). Interpretation:
**`rg` significantly below 1 ⇒ pervasive genotype × haplogroup interaction.** (`rg`≈1 with a wide SE
⇒ underpowered, not "no interaction"; below ~5 000 per stratum it's usually uninformative.)

**`h2_{I,R}.log`** — per-stratum SNP-heritability. Read the line
`Total Liability scale h2: <estimate> (<SE>)` (liability scale, using the population + per-stratum
sample prevalence). `Intercept` near 1 is good; `Intercept` ≫ 1 (and high `Lambda GC` / `Mean
Chi^2`) flags residual confounding/structure.

**`h2_full.log`** — the **pooled** (non-stratified, all-male) SNP-heritability, from a full-sample
main-effect GWAS (`gwas_full_*`). Same fields as `h2_{I,R}.log`.

**`h2_by_stratification.tsv`** *(when `H2_POOLED`)* — the **with- vs without-Y-haplogroup**
comparison in one table: one row each for `pooled`, `I`, `R`, with `h2`, `h2_se`, `intercept`,
`lambda_gc`, `mean_chi2`, `n_snps`. Compare the pooled h2 to the two stratum h2's; combined with a
cross-stratum `rg` < 1, a pooled estimate that differs from the stratum estimates points to
Y-haplogroup-dependent architecture. (The stratum SEs are wide on half-samples — read the
`h2_se`.)

**`h2_by_chromosome.tsv`** *(when `H2_PER_CHR`)* — per-chromosome SNP-heritability from the pooled
GWAS: one row per autosome `1`…`22` (and `X` when an X LD reference is supplied via `EUR_LD_X` +
`HM3_X`), so you can weigh each chromosome's contribution — e.g. **chrX vs the autosomes**. Each
row has the same columns as above; to compare contributions fairly look at h2 relative to each
chromosome's `n_snps` (per-chromosome h2 is noisy — small chromosomes especially — and these are
separate LDSC fits, the simple/interpretable choice; partitioned LDSC or GREML are more rigorous).

**`gwas_{I,R}_autism.regenie`** — the raw per-stratum main-effect GWAS (REGENIE step 2). Feeds Arm
B; `gwas_{I,R}_chr<c>_autism.regenie` chunks (if split) are intermediate.

**`gwas_{I,R}.forldsc.txt`** — those sumstats reshaped to LDSC-munge-ready columns
(`SNP A1 A2 N BETA SE P`; `LOG10P→P`, `ALLELE1→A1`). **`gwas_X_{I,R}.forldsc.txt`** are the chrX
equivalents — emitted ready to munge, but LDSC h2/rg on chrX needs an X-specific LD reference you
supply (not wired here).

**`munged_{I,R}.sumstats.gz`** — the LDSC-internal munged sumstats (`SNP A1 A2 Z N`) consumed by
`h2`/`rg`. Not meant for direct reading.

---

## LAVA — local cross-stratum r_g

LAVA runs in **two stages** (`METHODS.md` §7): a **genome-wide screen** (R, the headline)
that *nominates* blocks, then a **per-block deconfounding** (Python) for the blocks you
follow up. The screen is an effect-size landscape blind to ancestry; the deconfounding is
what makes a hit defensible.

### Genome-wide screen — `lava_local_rg.tsv` (when `LAVA_PARTITION` is set)

The headline local `r_g` from the real **LAVA R package**, **one row per LD block** in the
partition file (`blocks_s2500_m25_f1_w200.GRCh37_hg19.locfile`, ~2,500 blocks). Built by
`lava_inputs` → `lava_local`; needs R + the LAVA package (point `RSCRIPT` at it). Columns are
LAVA's bivariate output plus the block coordinates:

| column | meaning |
|--------|---------|
| `locus`, `chr`, `start`, `stop` | the LD block (coordinates in the genotype build, hg19) |
| `phen1`, `phen2` | the two strata correlated (`I`, `R`) |
| **`rho`** | **local cross-stratum genetic correlation** for the block — the headline; **well below 1, or sign-flipped, ⇒ the genetic architecture diverges between haplogroups here** |
| `rho.lower`, `rho.upper` | 95% CI for `rho` — a divergence candidate has `rho.upper` clearly **< 1** |
| `r2`, `p` | local explained correlation, and LAVA's p for `rho ≠ 0` |

> **Read `rho`, not `p`, for divergence.** LAVA's `p` tests `rho ≠ 0` (is there *any* local
> correlation); the I-vs-R question is whether `rho` is *below 1*. Nominate blocks by low
> `rho` / `rho.upper < 1`. Blocks with too few SNPs or no univariate signal are **silently
> dropped** by LAVA, so absence is not `rho = 1`. This screen does **not** deconfound —
> every nominee must go through the per-block stage below before you believe it.

### Per-block deconfounding — `lava_<name>_*` (when `LAVA_LOCI` is set)

One set of files per block you listed in `LAVA_LOCI` (the screen's nominees), named by its
label. This is the ancestry-matched permutation + negative-control panel.

**`lava_<name>_summary.txt`** — the headline for that block (`key<TAB>value`):

| key | meaning |
|-----|---------|
| `target_local_rg` | local cross-stratum genetic correlation at the block (sign matters) |
| **`target_anc_matched_p`** | ancestry-matched empirical p (small ⇒ survives the null) |
| **`target_vs_controls_tail_frac`** | where the block sits among SNP-count-matched **negative-control** blocks — **small ⇒ it stands out** from the genome-wide background |
| `control_local_rg_mean`, `control_p_median` | the control panel's centre |
| `strata_mean_size`, `within_stratum_HapPC_AUC` | the strata used |

A block is interesting only when **both** its `anc_matched_p` is small **and** its
`tail_frac` is small — otherwise residual structure depresses local `rg` everywhere (see
`METHODS.md` §7). Same global-PC ceiling as Arm A. **Pair it with the screen:** the
`lava_local_rg.tsv` `rho` is the effect size, this `anc_matched_p` is the deconfounded
significance.

**`lava_<name>_loci.tsv`** — the target row (`is_target=True`) plus every negative-control block,
each with `local_rg` and `anc_matched_emp_p`, so you can see the distribution behind the tail
fraction.

**`lava_<name>_strata_selection.tsv`** — the strata tightening curve, as in Arm A.

---

## Female negative control (when `FBFILE` is set)

Females carry no Y, so a real Y-haplogroup interaction can't exist in them — these files
re-test the male hits in an ancestry-matched female split to expose ancestry artefacts.

**`female_negative_control.tsv`** — one row per male interaction hit:

| column | meaning |
|--------|---------|
| `ID` | the variant (a male interaction hit) |
| `male_LOG10P` | its interaction `LOG10P` in the male scan |
| `female_int_chi2`, `female_int_p` | the `SNP × pseudo-Hap` interaction in the **females** |
| `looks_like_ancestry_artifact` | `True` if `female_int_p < α` (default 0.05) |

Read it asymmetrically: **`looks_like_ancestry_artifact = True` is strong evidence the male
hit is ancestry** (something with no Y shows the "interaction"). A **non-flagged** hit (large
`female_int_p`) is *consistent* with being Y-driven but not proof — autism's architecture is
sex-differential, so a hit can be null in females for sex reasons rather than because it needs
the Y. Treat positives as confounding confirmed, nulls as "not refuted." (See `METHODS.md` §6.)

**`female_lambda.txt`** — `females_used`, `hits_tested`, `hits_flagged_artifact`, and the
genome-wide `female_interaction_lambda` over a random panel (overall inflation of the female
SNP×pseudo-Hap scan). `female_pseudohap.txt` is the female pseudo-haplogroup assignment
(intermediate).

## A 30-second decision guide

- **Is a specific variant a real I-vs-R interaction?** `top_interactions(_X).tsv` genome-wide
  significant **and** `perm_interactions.tsv` `anc_matched_emp_p` small. Even then: provisional
  (global-PC ceiling) → confirm with local ancestry / replication.
- **Does haplogroup reshape the architecture overall?** `I_vs_R_rg.log`: `rg` significantly < 1.
- **Does a region diverge locally?** Two steps: screen `lava_local_rg.tsv` for blocks with `rho`
  well below 1 (`rho.upper < 1`), then deconfound each nominee — `lava_<name>_summary.txt` with
  small `target_anc_matched_p` **and** small `target_vs_controls_tail_frac`.
- **Is a surviving interaction hit actually ancestry?** `female_negative_control.tsv`: if it's
  flagged (`looks_like_ancestry_artifact`), yes — kill it.

## Everything else (not results)

`*.log` (qc/step1/regenie/ldsc/plink), `*_pred.list`, `*.loco`, `keep_{I,R}.txt`,
`int_covars.txt`, `qc_pass.snplist`, `lava_pcs.*`, `lava_prune.*`, the pooled-GWAS
intermediates (`gwas_full_*.regenie`, `gwas_full.forldsc.txt`, `munged_full.sumstats.gz`,
the per-chromosome `h2_chr<c>.log`), and the `*_chr<c>_*` / `*_b<k>_*` chunk/batch files are
**intermediate** — inputs to later steps or scratch for debugging, not final results.
