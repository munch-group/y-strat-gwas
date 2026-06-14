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
   - A **female cohort** (no Y, same SNP IDs as the males) whose phenotype keeps the confound
     SNP's g×ancestry effect but has *no* Hap interaction — so the negative control should flag
     the confound (reproduces in females) but not the real interaction.
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
   - per-stratum **chrX** summary stats (`gwas_X_{I,R}.forldsc.txt`) are emitted,
   - the **females negative control** flags the confound SNP (reproduces in a no-Y group) but not
     the genuine interaction SNP,
   - the **pooled vs stratified** and **per-chromosome** h2 tables are produced.

The `lava_local` headline (real LAVA R package) is **not** run by the test — R/LAVA isn't in the
pixi env. Only the pure-Python deconfounding layer (`lava_perm_*`) is exercised.

## Parallelisation checks (`test_parallel.py`)

After the main test (it reuses the data + intermediate results), run:

```bash
pixi run --manifest-path ./pixi.toml python tests/test_parallel.py
```

It verifies the two parallelisation features: chromosome-split REGENIE step 2 + `concat_regenie.py`
reproduces the unsplit scan **exactly**, and pooling K permutation batches with `pool_perm.py`
reproduces the single-job p-values (and sums the permutation count) for both `perm_interaction`
and the LAVA local-r_g test.

The synthetic genotypes span **chr 1–22** (~182 variants each), so the data splits across tasks
out of the box. `tests/env.sh` therefore sets `YS_SPLIT_CHROMS=1-22` and `YS_PERM_BATCHES=4`, so a
gwf run against the test data (`source tests/env.sh; gwf -b local run`, or `tests/run_via_gwf.sh`)
fans out into ~104 tasks: one step-2 job per chromosome for `interaction`/`gwas_{I,R}` (+ gathers)
and four pooled permutation batches each. The whole split pipeline has been verified to complete
through the real gwf daemon with correct gathered/pooled results. Drop those two env vars to run
the single-job path. (The in-process `run_pipeline_test.py` does not set them, so it stays
single-job and fast.)

Note on the gwf **local** backend: a fast-finishing wave of parallel jobs can leave the downstream
gather in `shouldrun` until the next `gwf run` (Slurm chains this via job dependencies);
`run_via_gwf.sh` re-submits automatically when the DAG stalls.

## Real-data input adapters (`test_realformat.py`)

The generator also emits `haplogroup_major.txt` (haplogroup in a column named `Major`, with ~15%
non-I/R males), `raw_covariates.txt` (`FID IID SOL C1..C10 st1`, with a header) and the
**headerless** `raw_phenotypes.txt` (`FID IID value`, 1/2-coded) — mirroring the real iPSYCH
file layouts. After the main test, run:

```bash
pixi run --manifest-path ./pixi.toml python tests/test_realformat.py
```

It checks `make_strata --hap-col Major` drops the non-I/R males (and emits `keep_IR`), `prep_covar`
renames `C→PC` and adds dummy age/batch, `make_interaction_covars` restricts to I/R, `prep_pheno`
recodes the **headerless** PLINK 1/2 `.pheno` → REGENIE 0/1 (`--no-header --value-col 3`), and that
with the `RAW_*` env set the gwf graph orders the prep tasks **before** `step1_full` (the prepped
files become `PHENO`/`BASECOVAR`).

## Notes / environment quirks handled here

- Genotypes are generated as a VCF and imported with `plink2 --id-delim _` so the
  `FIDk_IIDk` sample names split back into matching FID/IID columns.
- `KMP_AFFINITY=disabled` is exported because MKL/OpenMP thread-affinity binding asserts
  in some sandboxed environments; it's unrelated to the pipeline.
- Estimates (h2, rg, λ_GC) are **noise at this scale** — the test proves the stages run,
  interoperate, and point the right way, not that the numbers are precise.
