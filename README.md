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

## Input formats

| file | columns |
|------|---------|
| `genotypes.{bed,bim,fam}` | plink1, autosomes, QC'd samples (X handled separately, see below) |
| `phenotypes.txt` | `FID IID autism` — `autism` = 1 case / 0 control / `NA` |
| `covariates_base.txt` | `FID IID PC1 … PCk age batch …` |
| `haplogroup.txt` | `FID IID Hap` — `Hap` = `I` or `R` |

`Hap` is coded internally as I=1, R=0.

## Setup

```bash
# modern tools
pixi install                          # reads pixi.toml (regenie, plink2, py3, gwf)

# classic LDSC (separate python 2.7 env)
mamba env create -f environment-ldsc.yml
git clone https://github.com/bulik/ldsc.git    # -> set LDSC_DIR in workflow.py
```

Download LDSC reference data (distributed by the LDSC authors; mirrored on
Zenodo) and set the paths in `workflow.py`:

- `eur_w_ld_chr/` — 1000G EUR LD scores → `EUR_LD`
- `w_hm3.snplist` — HapMap3 SNPs for `--merge-alleles` → `HM3`

## Configure & run

Edit the CONFIG block at the top of `workflow.py` (paths, `ACCOUNT`, `NPC`,
`CATCOVAR`, `PREV_POP`, `STRATUM_SPECIFIC_STEP1`).

```bash
gwf -f workflow.py status
gwf -f workflow.py run
```

Target graph: `qc`, `strata`, `int_covars` → `step1_full` →
`interaction` → `top_int`; and per stratum `gwas_{I,R}` → `munge_{I,R}` →
`h2_{I,R}`, `rg`.

## Reading the results

- `results/top_interactions.tsv` — strongest SNP×Hap signals, plus `lambda_GC`
  on the interaction p-values printed to the log. If λ > ~1.10, revisit QC and
  the PC×Hap covars before trusting hits. Genome-wide interaction threshold is
  5e-8 unless you used the two-step screen (below).
- `results/I_vs_R_rg.log` — the cross-stratum r_g (parsed to stdout by `rg`).
- `results/h2_{I,R}.log` — per-stratum SNP-h2.

## Important caveats baked into / around this pipeline

1. **Structure is the main threat.** PC×Hap terms guard the E-side; REGENIE
   does not add per-SNP SNP×PC (G-side) terms, so for any genome-wide-significant
   interaction hit, re-test it within an ancestry-homogeneous subset or with
   local-ancestry covariates before believing it.
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
5. **X chromosome.** Excluded here (LDSC is autosomal; male hemizygosity needs
   special coding). Given your X-linked interest, run X as a separate REGENIE
   step-2 pass with appropriate male coding.
6. **Mechanism.** The binary I/R label is a coarse proxy. A continuous
   Y-heterochromatin estimate (relative read depth over DYZ1/DYZ2) may capture
   the sink effect better — swap it in as the interacting covariate by replacing
   `Hap` with that column.
