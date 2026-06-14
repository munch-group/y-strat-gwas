# Detecting SNP × Y-haplogroup epistasis under ancestry confounding

*A methodological and practical guide to the analysis this repository implements.*

This document is meant to be read on its own. It introduces the statistical ideas
needed to test for interaction between Y-chromosome haplogroup and the rest of the
genome, explains the one problem that dominates everything (haplogroup is not an
exogenous variable — it tags ancestry), and shows the concrete approaches the
pipeline uses to confront it. The companion `README.md` documents how to *run* the
workflow; this document explains *why it is built the way it is*.

---

## 1. The question

In a cohort of male autism cases and controls, the Y chromosome splits the men into
two deep lineages — haplogroup **I** and haplogroup **R**. The scientific question is
whether the effect of autosomal (and X-linked) genetic variation on autism *depends
on* which Y-haplogroup a man carries. Two flavours of the question:

- **Individual variants.** Are there specific loci whose association with autism
  differs between I and R men? (A genome-wide SNP × haplogroup interaction scan.)
- **General architecture.** Does haplogroup reshape the genetic architecture as a
  whole — i.e. is there *pervasive*, polygenic genotype × haplogroup interaction,
  even if no single variant reaches significance? This is asked at two resolutions:
  one genome-wide number (cross-stratum genetic correlation) and one per LD block
  (local genetic correlation, for following up a specific region).

These map onto the pipeline's analysis arms below — the individual-variant scan
(Arm A, Section 3), the genome-wide architecture readout (Arm B, Section 4) and its
local-resolution companion (LAVA, Section 7), with the X chromosome handled
separately (Section 8). Each rests on a shared deconfounding idea (Section 6), and
each is easy to state and treacherous to estimate, for one reason.

---

## 2. The dominant problem: haplogroup is not exogenous

A clean interaction analysis wants the "environment" — here, haplogroup — to be
*independent* of the rest of the genome except through the biology you care about. Y
haplogroup badly violates this. In Northern Europe, the Y lineages co-segregate with
**autosomal ancestry**: R (R1b/R1a) men carry more Steppe/Yamnaya-related ancestry,
while I (I1/I2) men skew toward earlier Western-Hunter-Gatherer / Early-European-
Farmer components. So "I vs R" is, to a substantial degree, a **proxy for a
continuous autosomal ancestry gradient**.

This single fact reorganises the whole analysis. It means that **any autosomal
variant whose frequency or effect differs along that ancestry gradient can produce a
SNP × haplogroup signal that has nothing to do with epistasis.** A balanced 50/50
split of the cohort does not help; the confounding is in the *correlation structure*,
not the marginal counts.

### How a spurious interaction arises

It is worth being precise about the mechanism, because the fix follows directly from
it. A fake SNP × Hap interaction is generated whenever:

1. a SNP's **effect on the trait varies along ancestry** (a genuine SNP × ancestry
   interaction), **and**
2. **Hap is correlated with ancestry** (which it is, by construction).

Crucially, you do **not** need any real biological epistasis for (1) — it arises
mechanically from **differential linkage disequilibrium**. A genotyped tag SNP is in
different LD with the true causal variant in a Steppe-rich versus a WHG-rich
background, so its *estimated* per-allele effect differs across the gradient. Because
Hap tags that gradient, a model that can only attribute effect-modification to Hap
will load this SNP × ancestry signal onto the SNP × Hap term. The result is a
confident, entirely artefactual "interaction."

### Global versus local ancestry

There are two resolutions of ancestry, and the distinction is central later:

- **Global ancestry** — a person's genome-wide ancestral composition, summarised by
  the top principal components (PCs) of the genotype matrix.
- **Local ancestry** — the ancestral origin of the specific chromosomal segment *at a
  given locus*, which fluctuates along the genome and is only imperfectly summarised
  by global PCs.

The confounder that corrupts a *specific* hit is **local** ancestry at that locus.
Everything built on global PCs — covariate adjustment, PC-matching, the permutation
test below — controls confounding only to global-PC resolution and cannot, even in
principle, remove a purely local-ancestry artefact. Hold onto this; it is the
honest limit of the entire analysis.

---

## 3. Arm A — the individual-variant interaction scan

### Why a mixed model

Relatedness and population structure inflate naive association tests. The standard
modern solution is a **whole-genome regression / mixed-model** approach (the pipeline
uses REGENIE): a first step fits a polygenic predictor of the phenotype from a pruned
set of common variants (absorbing relatedness and structure into a per-person
offset), and a second step tests each variant conditional on that offset. For binary
traits with case/control imbalance, **Firth** or **saddle-point (SPA)** corrections
keep the test calibrated for rare variants.

### The interaction model, and the confounder terms

For a tested SNP `g` and haplogroup `Hap`, the model is

```
autism ~ g + Hap + g:Hap + PCs + Hap:PCs (+ covariates)
```

and the test of interest is the coefficient on `g:Hap`.

The non-obvious requirement is the `Hap:PCs` block. This is **Keller's (2014) point**
about confounded gene × environment studies: when the "environment" is correlated
with a covariate, it is not enough to include the covariate as a main effect — you
must also include the covariate's *interactions*. Concretely, to test `g × E`
without bias when a covariate `C` confounds it, the model needs both `C × E` and
`C × g` terms:

- **`C × E` = `PC × Hap`** guards against ancestry modifying the **Hap** effect.
- **`C × g` = `PC × g` (SNP × PC)** guards against ancestry modifying the **SNP**
  effect — which, as we just saw, is *exactly* the channel that manufactures spurious
  SNP × Hap signals.

### Why only one of the two guards is in the scan

There is a structural asymmetry between the two guard terms:

- **`PC × Hap` is per-person constant.** For a given individual it does not depend on
  which SNP is being tested, so all 20 columns can be **precomputed once** and dropped
  into the covariate file. The pipeline does exactly this
  (`make_interaction_covars.py`), and the interaction routine carries them as ordinary
  additive covariates. The E-side guard is therefore essentially free.
- **`SNP × PC` is per-SNP.** It changes for every variant in the scan, so it can only
  be fit *inside* the per-variant model. REGENIE's interaction machinery tests a
  *single* interaction term (`g × Hap`) and crosses `g` with nothing else; it does not
  — and is not built to — add `g × PC1..PCk` for every variant. So the G-side guard is
  simply absent from the genome-wide pass.

This is not an oversight you can flip a flag to fix; it is what the optimised
single-interaction test does. The consequence: **a genome-wide-significant SNP × Hap
hit is provisional until the missing G-side guard is supplied some other way** — which
is the job of Arm A's post-hoc permutation re-test (Section 6).

### Power is the binding constraint

Detecting an interaction of a given per-allele magnitude needs roughly **four times**
the sample size of detecting a main effect of the same magnitude (the interaction
estimate has larger variance). A naive genome-wide interaction scan at 5×10⁻⁸ will
usually find nothing. Three responses, in increasing order of how much they help:

- **Two-step screening.** First filter variants on a quantity that is independent of
  the interaction under the null — marginal variance heterogeneity (a vQTL pre-scan),
  or marginal G–Hap association (Murcray's approach) — then test the interaction only
  in the survivors, slashing the multiple-testing burden.
- **Candidate-set restriction.** Up-weight or restrict to loci where a mechanistic
  prior predicts signal (here: imprinted clusters, X-linked regulators). This both
  cuts the testing burden and, importantly, gives a prior that is **orthogonal to the
  ancestry-artefact process** — a hit in a pre-specified mechanistic set is more
  credible than a random genome-wide hit.
- **Do not use the case-only design.** Testing G⊥Hap within cases is tempting because
  it is far more powerful, but it *assumes* G and Hap are independent in the
  population. Ancestry makes that assumption false here, so the case-only test is
  invalid in exactly this setting.

---

## 4. Arm B — the well-powered global readout

Single-variant interaction tests are underpowered; the *aggregate* question is much
better powered. Treat **"autism among I men"** and **"autism among R men"** as two
separate traits and ask how genetically similar they are.

### Cross-stratum genetic correlation

Run a main-effect GWAS *within* each stratum, then estimate the **genetic
correlation** `r_g` between the two sets of summary statistics. If haplogroup did not
reshape the architecture at all, the same alleles would have the same effects in both
strata and `r_g = 1`. **`r_g` significantly below 1 is direct evidence of pervasive
genotype × haplogroup interaction** — and because it pools information across the whole
genome, it is far better powered than any single-variant test.

### LD Score Regression, briefly

The pipeline estimates `r_g` (and per-stratum SNP-heritability) with **LD Score
Regression (LDSC)**. The idea: under polygenicity, a variant's association
test-statistic scales with how much genetic variation it *tags*, measured by its **LD
score** (the sum of squared correlations with neighbouring variants). Regressing
test-statistics on LD scores separates true polygenic signal (the slope, ∝
heritability) from inflation due to confounding or structure (the intercept).
The bivariate version yields the genetic covariance, hence `r_g`. LDSC needs only
summary statistics plus a reference panel of LD scores, which makes it cheap at
biobank scale.

Two practical points the pipeline handles:

- **Liability scale.** For a case/control trait, SNP-heritability must be converted to
  the underlying liability scale, which requires **both** the population prevalence and
  the **sample** (case) prevalence. LDSC errors if given only one; the pipeline
  computes the per-stratum case fraction (`samp_prev.py`) and supplies both.
- **Shared vs stratum-specific step 1.** Reusing the full-sample polygenic predictor
  for both per-stratum GWASs is cheaper but biases `r_g` slightly **toward 1** (i.e.
  conservative for the hypothesis). Fitting step 1 separately within each stratum is
  cleaner at ~2× the compute. Either way, with the cohort halved the LDSC standard
  errors are the binding constraint — below roughly 5 000 per stratum the `r_g` is
  unlikely to be informative, and that is the first number to sanity-check.

LDSC as distributed is autosomal, so Arm B's `r_g`/`h2` do **not** include the X
chromosome (see Section 8).

### Heritability: stratified vs pooled, and per-chromosome

Two complementary readouts sit alongside `r_g`, both from LDSC `--h2`.

**Stratified vs pooled (does ignoring haplogroup change the apparent heritability?).**
The pipeline reports per-stratum SNP-heritability `h2_I`, `h2_R` and a **pooled**
estimate `h2_pooled` from a single full-sample GWAS that ignores Hap (all males
together). If haplogroup does not reshape the architecture, the three agree (the
pooled estimate is the tightest, since it uses the whole sample). If there is genuine
G × Hap interaction, the alleles' effects differ between strata, so the pooled GWAS
estimates an *average* effect — and when effects are heterogeneous they partially
cancel, which can pull `h2_pooled` **below** the stratum estimates. Read this as a
**descriptive companion to `r_g`**, not an independent test: `r_g < 1` is the powered,
rigorous statement; the h2 comparison is a sanity check in the same direction. And
mind the SEs — each stratum is a half-sample, so its h2 SE is wide and apparent
pooled-vs-stratum gaps are often noise (`h2_by_stratification.tsv` carries the SE for
exactly this reason).

**Per-chromosome (how much does each chromosome — especially chrX — contribute?).**
SNP-heritability can be partitioned by chromosome by running LDSC once per chromosome
(restricting both the sumstats and the LD scores to that chromosome). Under
polygenicity each chromosome's h2 is roughly proportional to its length / SNP count,
so the interesting question is **enrichment**: does a chromosome carry more (or less)
heritability than its size predicts? Given this project's interest in X-linked
regulators, the headline use is **chrX versus the autosomes** — compare h2 *per SNP*
(h2 ÷ `n_snps`), not raw h2, to put a small chromosome on equal footing. Caveats:
per-chromosome estimates are **noisy** (few SNPs, especially small chromosomes), these
are separate fits rather than one joint model (partitioned LDSC or GREML with
per-chromosome GRMs are more rigorous), and — because LDSC is autosomal — **chrX is
included only when you supply an X-chromosome LD reference** (the pooled chrX GWAS is
emitted ready to feed it; see Section 8).

---

## 5. Interaction is scale-dependent

A statistical interaction can appear on one scale and vanish on another. A SNP × Hap
effect that is significant on the **logit** scale (the natural scale of the
logistic model) may disappear on the **liability** scale, and vice versa. This is not
a nuisance to be smoothed over — it is a reason to weight **qualitative** interactions
(where the *direction* of an effect flips between I and R) above **quantitative** ones
(where only the magnitude differs and can often be removed by a monotone
re-scaling). Before attaching any biological interpretation to an interaction, check
that it is robust across scales.

---

## 6. The post-hoc G-side guard: ancestry-matched permutation

Section 3 left a gap: the scan cannot include the per-SNP `SNP × PC` term, so a hit
might be a real I-vs-R effect or an ancestry artefact, and the two are
indistinguishable from the scan output alone (they can have near-identical test
statistics). The pipeline closes this gap for the hits with an **ancestry-matched
permutation test** (`ancestry_matched_perm.py`, target `perm_interaction`). This is
the one piece of new methodology worth understanding in detail.

### The idea

We want a null distribution for the SNP × Hap statistic that **preserves the
Hap–ancestry confounding but removes any genuine epistasis**, so that beating it means
"more than ancestry alone can explain." We get it by **permuting Hap within
ancestry-matched strata**:

- Group individuals into strata that are **homogeneous in ancestry**.
- Within each stratum, **shuffle the Hap labels** among members, holding genotype,
  phenotype and covariates fixed.
- Recompute the SNP × Hap statistic on the shuffled labels; repeat many times to build
  the null; compare the observed statistic to it.

### Why this separates real interaction from confounding

The key is *where* each kind of signal lives relative to the strata.

- A **genuine** `g × Hap` effect is a **within-stratum** phenomenon: even among men of
  the same ancestry, the SNP's effect differs between I and R. Permuting Hap within
  the stratum destroys this association, so the observed statistic **exceeds** the
  permutation null ⇒ small empirical p ⇒ **survives**.
- A **pure ancestry artefact** is a **between-stratum** phenomenon: within a stratum
  ancestry is ~constant, so the SNP's effect is ~constant and genotype is independent
  of Hap; the artefact's signal lives entirely in the differences *between* strata,
  which within-stratum permutation **preserves**. The observed statistic therefore
  sits **inside** the null ⇒ empirical p → 1 ⇒ **collapses**.

So the permutation acts as a sieve: it leaves within-ancestry interaction intact in
the contrast between observed and null, while folding the between-ancestry confound
into the null itself.

### How to build the strata: propensity, not raw PC distance

The strata are only as good as their ancestry resolution, and there is a subtle trap.
The naive choice — cluster individuals by Euclidean distance in PC space — **fails**
when most PCs are ancestry-uninformative: the clustering wastes its resolution on
noise dimensions and never actually separates the ancestry axis that matters. In a
realistic case you can drive the cluster sizes down and *still* find that Hap is
predictable from PCs within each cluster.

The principled fix is **propensity-score stratification** (Rosenbaum–Rubin). Fit the
**Hap propensity** `π = P(Hap = I | PCs)` — a single direction in PC space, the one
along which I and R actually separate — and bin individuals by `π`. Within a bin,
everyone has the same probability of being I given ancestry, so by construction Hap is
~independent of ancestry, and the bin is "ancestry-matched" in precisely the sense the
confound requires. Because the propensity is one-dimensional, binning concentrates all
the resolution on the confounding axis and ignores the noise PCs automatically.

### When are the strata "homogeneous enough"? A principled stopping rule

"Ancestry-homogeneous" does **not** have to mean "everyone genetically identical." For
breaking the confound it means the weaker, cheaper condition: **within a stratum, Hap
is no longer predictable from the PCs** (`Hap ⫫ ancestry | stratum`). That is directly
measurable and gives an operational definition rather than a vibe:

> Regress Hap on the within-stratum-centred PCs and compute the AUC. **AUC ≈ 0.5**
> means the PCs carry no residual information about Hap inside strata — the target.

This yields an automatic tightening rule. Sweep stratum granularity from loose to
fine; for each, compute the within-stratum Hap~PC AUC; then take the **loosest** strata
(largest, for maximal power and permutation entropy) that still reach AUC ≈ 0.5.
Tightening further only sheds power without reducing detectable confound; loosening
re-introduces it. The pipeline writes this whole selection curve to
`perm_strata_selection.tsv` so the operating point is auditable.

### What it can and cannot do

This is a faithful implementation of the matched-subset idea using the **whole**
sample (every stratum contributes, rather than carving out one homogeneous slice), with
an exact permutation null. But it inherits the hard limit from Section 2:

> It controls confounding **only to global-PC resolution.** A hit driven by **local
> ancestry** at its own locus — genotype differing between I and R beyond what global
> PCs capture — produces a within-stratum association that the permutation cannot
> remove, and it will survive.

So a survivor means **"not explained by global ancestry structure,"** not "proven
epistasis." The genuinely clean confirmations need information this cohort does not
contain — **local-ancestry calls** at the locus (RFMix/FLARE-style), or **replication**
in an independent sample. The permutation test is the strongest filter available from
the data in hand, and survivors should be reported as provisional pending one of those.

### A genome-wide companion diagnostic

The same machinery yields a scan-level check: compare the **observed** distribution of
interaction statistics (its genomic-inflation factor λ) to the **permutation-null** λ
(`perm_lambda.txt`). This is a far more honest calibration reference than comparing λ
to 1, because the null is built against the actual Hap–ancestry structure rather than
an idealised no-structure world.

### A no-Y negative control: females

There is a second, independent line of defence that uses a different population rather
than a different null: **females**. A female carries no Y chromosome, so a genuine
Y-haplogroup × autosome interaction *cannot* exist in her — but the **ancestry artefact
can**, because females sit on the same autosomal ancestry gradient that haplogroup I vs
R tags in males. So we split the females along exactly that axis and look for the
interaction there:

1. fit the Hap **propensity** `P(Hap = I | PCs)` on the males and apply it to the
   females, labelling each female pseudo-I / pseudo-R to match the male I:R ratio
   (`assign_pseudo_hap.py`) — this is the *same* ancestry split, transplanted into a
   group with no Y;
2. re-test each **male interaction hit** for `SNP × pseudo-Hap` in the females, with the
   same model (`female_negcontrol.py`).

The readout is sharp: **a male hit that reproduces in the females is an ancestry
artefact** (something with no Y is showing the "interaction"); a male hit that stays
**null** in the females is consistent with being genuinely Y-driven.

Two things make this complement the permutation rather than duplicate it. First, the
female split is not conditioned on PCs at the testing stage, so it can expose
artefacts driven by ancestry structure that the (global-PC-bounded) permutation test
misses — including *local* ancestry that tracks the global gradient. Second, and
crucially, the test is **asymmetric in what it proves**: a *significant* female result
is strong, direct evidence of confounding; a *null* female result is weaker — it is
consistent with a Y-driven effect, but autism's architecture is **sex-differential**,
so a hit could be null in females for reasons of sex rather than because it needs the Y.
Read a positive female result as "this hit is ancestry," and a null as "not refuted —
still needs replication / local ancestry." It also assumes the female PCs live in the
same ancestry space as the male PCs (PCs computed jointly, or females projected).

This is the principled way to bring females in: as a **negative control**, *not* as a
third symmetric stratum alongside I and R — pooling the sexes into one stratification
would let **sex** (a far larger axis than Y-haplogroup, with its own differential
architecture and X-dosage differences) dominate every contrast, reintroducing exactly
the kind of confounding the male-only design was built to avoid.

---

## 7. Localising the divergence: LAVA local genetic correlation

Arm B gives one genome-wide `r_g`; the single-variant scan gives individual loci.
Between them sits a third, valuable resolution: **local** genetic correlation within
each of ~2,500 LD blocks (the method here is **LAVA**). Treating ASD-in-I and
ASD-in-R as two traits, a block where the local `r_g` falls well below 1 — or flips
sign — is a region where the genetic architecture diverges between haplogroups, even
if no single variant in it reaches interaction significance. It is the natural way to
follow up a mechanistic candidate region (e.g. a specific gene).

**It is confounded by exactly the same problem, and arguably more visibly.** The
per-stratum effect estimates differ partly because they are measured in samples of
different ancestry composition, so **differential LD** depresses the cross-stratum
concordance with no epistasis required. A local `r_g` that sits below 1 *across the
whole genome* is precisely the fingerprint of residual structure. So a focused
local-`r_g` claim needs the same deconfounding as a single-variant hit — and one
addition.

The pipeline pairs the **headline LAVA estimate** (the canonical R tool, run on the
per-stratum summary statistics) with a **deconfounding layer** that mirrors Section 6
and extends it:

- **Same ancestry-matched permutation.** A local statistic — an LD-aware cross-stratum
  correlation of the per-stratum block betas, `T = bᵢ'R bᵣ / √(bᵢ'R bᵢ · bᵣ'R bᵣ)`,
  PC-residualised so global ancestry is regressed out — is recomputed under Hap
  permuted within Hap-propensity strata. A genuine within-ancestry divergence survives;
  an ancestry artefact, living between strata, collapses.
- **Region-excluded PCs (anti-circularity).** The ancestry features the locus is
  matched on are computed with the test locus's chromosome **excluded**, so the locus
  cannot leak into its own confounder adjustment.
- **A negative-control-loci panel.** The *identical* matched-permutation procedure is
  run at a panel of SNP-count-matched control blocks elsewhere in the genome. This is
  the decisive check, because a ratio statistic like local `r_g` can be depressed
  *genome-wide* by residual structure: the locus is only interesting if it sits in the
  **tail of the control distribution**, not merely if its single-locus p is small.
  Without the control panel you cannot distinguish "this locus is special" from "I/R
  local `r_g` is negative everywhere," and the latter is exactly what residual
  structure produces.

The same hard limit from Section 2 still applies — everything is built on global PCs,
so a **local-ancestry** artefact at the locus survives both the permutation and the
control panel (the controls narrow this to "not a genome-wide effect" but cannot rule
out a locus-specific local-ancestry artefact). A surviving, control-beating local `r_g`
is the strongest statement the data alone support; confirming it as biology still needs
local ancestry at the locus or replication.

## 8. The X chromosome

The X is mechanistically the most interesting target here (it hosts X-linked
regulators of interest), but it needs special handling for two reasons:

- **Hemizygosity.** Males carry one X, so non-pseudoautosomal (non-PAR) X genotypes
  must be coded hemizygously (0/2), using the sex of each individual and the build's
  PAR boundaries. REGENIE does this automatically when the chromosome is coded as X/23,
  the `.fam` sex column is filled, and `--par-region <build>` is given; the pipeline
  runs a dedicated chrX step-2 pass on a separate fileset (`XBFILE`), reusing the
  autosomal polygenic predictor.
- **LDSC is autosomal.** The bundled LD-score reference does not cover the X, so Arm
  B's `r_g`/`h2` exclude it. The pipeline still produces per-stratum chrX summary
  statistics in LDSC-ready form; running LDSC on them additionally requires an
  X-chromosome LD-score reference, which is left to the user.

The same ancestry-confounding logic, and the permutation guard, apply to chrX hits as
to autosomal ones.

---

## 9. How the concepts map to the pipeline

| Concept (this document) | Where it lives |
|---|---|
| Mixed-model interaction scan, Firth/SPA | `interaction` target (REGENIE step 2, `--interaction Hap`) |
| E-side guard `PC × Hap` (Keller) | `int_covars` target / `make_interaction_covars.py` |
| Missing G-side guard `SNP × PC` | *deliberately absent* from the scan — see next row |
| Ancestry-matched permutation re-test (G-side guard) | `perm_interaction` / `ancestry_matched_perm.py` |
| Propensity stratification + AUC tightening rule | `ancestry_matched_perm.py` → `perm_strata_selection.tsv` |
| Genome-wide λ vs permutation null | `perm_lambda.txt` |
| Cross-stratum genetic correlation `r_g` | `gwas_{I,R}` → `munge_{I,R}` → `rg` (LDSC) |
| Local cross-stratum `r_g` (headline) | `lava_inputs` → `lava_local` (LAVA, R) |
| Local-`r_g` deconfounding + negative controls | `lava_perm_<name>` / `local_rg_perm.py` |
| Region-excluded PCs (anti-circularity) | `lava_pcs` (PCA with the locus chromosome dropped) |
| Per-stratum liability-scale `h2` (samp+pop prev) | `h2_{I,R}` / `samp_prev.py` |
| Stratified vs pooled `h2` | `gwas_full` → `h2_full`, tabulated in `h2_by_stratification.tsv` (`collect_h2.py`) |
| Per-chromosome `h2` (chrX vs autosomes) | `h2_by_chromosome` / `collect_h2.py` (chrX needs `EUR_LD_X`+`HM3_X`) |
| Shared vs stratum-specific step 1 | `STRATUM_SPECIFIC_STEP1` |
| chrX hemizygous coding | `interaction_X` / `gwas_X_{I,R}` (`--par-region`, `--sex-specific male`) |

---

## 10. Limitations, stated plainly

- **The ceiling is global-PC resolution.** PC adjustment, PC×Hap terms, PC-propensity
  matching and the permutation test all control ancestry confounding only as far as
  the global PCs see it. **Local ancestry at a locus is never removed.** Treat every
  surviving hit as provisional pending local ancestry or replication.
- **Power.** Halving the cohort by haplogroup and the intrinsic ~4× interaction
  penalty mean both arms are SE-limited. Below ~5 000 per stratum the `r_g`/`h2` are
  unlikely to be informative; the genome-wide single-variant scan needs either the
  two-step screen or a candidate-set restriction to have a chance.
- **Haplogroup is a coarse proxy.** The binary I/R label is a noisy stand-in for the
  mechanistically relevant variable (e.g. a continuous Y-heterochromatin dosage). A
  continuous measure could be swapped in wherever `Hap` appears for a more direct test.
- **Scale dependence.** Prefer sign-flipping (qualitative) interactions and verify
  robustness across the logit and liability scales before claiming biology.

---

### Key references (entry points, not exhaustive)

- Keller MC (2014), *Gene × Environment interaction studies have not properly
  controlled for potential confounders*, **Biol Psychiatry** — the `C×E`/`C×G` rule.
- Murcray CE et al. (2009) — two-step gene–environment screening.
- Bulik-Sullivan BK et al. (2015), LD Score Regression (h²) and cross-trait LDSC
  (`r_g`) — **Nat Genet**.
- Werme J et al. (2022), LAVA — local genetic correlation — **Nat Genet**.
- Rosenbaum PR & Rubin DB (1983) — propensity-score stratification.
- Mbatchou J et al. (2021), REGENIE — **Nat Genet**.
