# Local end-to-end test

A self-contained smoke test that runs the **entire** GxHaplogroup workflow on tiny
synthetic data — no Slurm, no downloads, no real cohort.

```bash
pixi run --manifest-path ./pixi.toml python tests/run_pipeline_test.py
```

Runtime is a few minutes on a laptop. Everything lands in `tests/work/` (git-ignored);
delete it to start clean.

## What it does

1. **`make_dummy_data.py`** generates a synthetic dataset into `tests/work/data/`:
   - 2000 males (1000 haplogroup **I**, 1000 **R**), 4000 SNPs across chr 1–22.
   - Haplogroup is deliberately **confounded with ancestry** (an ancestry axis read out
     by PC1 differs by group, and allele frequencies drift along it) — the exact problem
     the `PC×Hap` covariates in Arm A exist to absorb.
   - **One planted autosomal and one planted chrX SNP×Hap interaction** (effect present in I,
     absent in R), plus **one confound-only SNP** (effect varies with ancestry, no true Hap
     interaction) — all recorded in `truth.txt`.
   - **Two LD blocks for the LAVA arm**: a *genuine* local divergence (chr20, effects flip sign
     by Hap) and a *confound* one (chr21, effects modulated by ancestry) — also in `truth.txt`.
   - A separate **chrX fileset** (`genoX.vcf` + `sex.txt`): all-male, hemizygous genotypes
     (dosage 0/2) at non-PAR positions, so REGENIE applies male hemizygous coding.
   - Matching `phenotypes.txt`, `covariates_base.txt`, `haplogroup.txt`, a HapMap3-style
     `w_hm3.snplist`, and a synthetic LDSC reference (`eur_w_ld_chr/`, real LD scores
     computed from the dummy genotypes — so the test never needs LDSC's bed reader).

2. **`run_pipeline_test.py`** imports `workflow.py` with the inputs supplied via `YS_*`
   env vars and executes **each gwf target's actual `spec`** in dependency order (the same
   shell commands gwf would submit), asserting every declared output appears.

3. Validation checks the result is not just well-formed but *correct in direction*:
   - the planted autosomal interaction SNP lands in the **top 5%** of the SNP×Hap scan,
   - the planted **chrX** interaction SNP lands in the **top 10%** of the hemizygous chrX scan,
   - the **ancestry-matched permutation** separates the real interaction (survives, low
     `anc_matched_emp_p`) from the confound-only SNP (collapses), and the strata decouple Hap
     from PCs (within-stratum AUC ≈ 0.5) — the two often have *near-identical* naive interaction
     statistics, so this is the check that the deconfounding actually works,
   - the cross-stratum **rg** summary is produced,
   - per-stratum **liability-scale h2** is reported,
   - the **LAVA** local-r_g permutation separates the genuine local divergence (survives, low
     `anc_matched_emp_p`, in the tail of the negative-control loci) from the confound block
     (collapses) — the deconfounding check for the local-r_g arm,
   - per-stratum **chrX** summary stats (`gwas_X_{I,R}.forldsc.txt`) are emitted.

The `lava_local` headline (real LAVA R package) is **not** run by the test — R/LAVA isn't in the
pixi env. Only the pure-Python deconfounding layer (`lava_perm_*`) is exercised.

## Notes / environment quirks handled here

- Genotypes are generated as a VCF and imported with `plink2 --id-delim _` so the
  `FIDk_IIDk` sample names split back into matching FID/IID columns.
- `KMP_AFFINITY=disabled` is exported because MKL/OpenMP thread-affinity binding asserts
  in some sandboxed environments; it's unrelated to the pipeline.
- Estimates (h2, rg, λ_GC) are **noise at this scale** — the test proves the stages run,
  interoperate, and point the right way, not that the numbers are precise.
