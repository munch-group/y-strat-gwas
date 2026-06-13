In my Danish autism cohorts, males split evenly into two groups carrying a Y chromosome from either haplogroup I or R. I want to investigate if there are any epistatic interactions between Y chromosome haplogroup and variants in the remaining genome. I am interested in individual variants and in general patterns og heritability related to haplogroup


The thing that will dominate this analysis, before any method choice, is that haplogroup I vs R is not a clean exogenous variable — it tags deep autosomal ancestry. In Denmark, R (R1b/R1a) carries more Steppe/Yamnaya-related autosomal ancestry while I (especially I1/I2) skews toward earlier WHG/EEF components. So any autosomal variant whose frequency differs across those ancestry gradients will produce a spurious SNP×haplogroup signal, and a halved-and-balanced design does nothing to remove it. Everything below is really about separating genuine epistasis from this structure.

**Individual-variant scan**

Use a mixed-model interaction GWAS (REGENIE's `--interaction` for binary traits, with Firth/SPA; SAIGE-GENE if you also want set-based) rather than plain logistic regression, so relatedness and structure are absorbed. The non-obvious requirement is Keller's (2014) point: when the "environment" is confounded, you must put the confounder's interactions in the model, not just its main effect. So fit

`autism ~ SNP + Hap + SNP:Hap + PCs + Hap:PCs (+ SNP:PCs) + covariates`

and test `SNP:Hap`. Omitting `Hap:PCs` lets ancestry leak directly into the interaction estimate. The `Hap:PCs` block is per-person-constant, so it's cheap; only `SNP` and `SNP:Hap` (and optionally `SNP:PCs`) are per-SNP.

Power is the real constraint — an interaction needs roughly 4× the N of an equivalent main effect. Don't run a naive genome-wide interaction scan at 5e-8 and expect anything. Use a two-step filter: screen on a quantity independent of the interaction under the null (marginal variance heterogeneity, or marginal G–Hap association à la Murcray, or a vQTL pre-scan), then test the interaction only in the survivors with a far lighter multiple-testing burden. Avoid the case-only design here even though it's tempting for power — it assumes G ⊥ Hap in the population, which your ancestry confounding violates outright.

**General pattern / heritability**

The cleanest global test of "does haplogroup reshape the architecture" is the cross-stratum genetic correlation. Treat autism-in-I and autism-in-R as two traits and estimate r_g, either by bivariate GREML or, more robustly at iPSYCH scale, by running a GWAS within each stratum and feeding the two summary stats to cross-trait LDSC. r_g significantly below 1 is direct evidence of pervasive haplogroup × autosome interaction, and it's far better powered than any single-variant test. Complement it with GCTA `--gxe` (V_G, V_GxHap, V_e) and with stratified h²_SNP estimates per group. For all of these on case/control, convert to the liability scale using the male-specific population prevalence, since the cohort is all male.

**Aligning with your heterochromatin-sink prior**

Two things buy you power and biological interpretability. First, the binary I/R label is probably a noisy proxy for the mechanistically relevant variable — Y heterochromatin dosage. Estimate a continuous Y-heterochromatin measure from sequencing depth over the DYZ1/DYZ2 (Yq12) repeats and use that as the interacting covariate; it may capture a sink effect that the discrete haplogroup smears out. Second, rather than weighting all autosomal variants equally, restrict or up-weight the interaction scan to loci your model predicts are heterochromatin-sensitive — imprinted clusters and the X (hemizygous in males, and the home of your X-linked MAPT regulators like HDAC6). A candidate-set interaction test over imprinted/X/your meiosis gene sets has dramatically better power than genome-wide and is the natural first pass given your hypothesis.

One caveat to keep in mind throughout: statistical interaction is scale-dependent. A signal on the logit scale can vanish on the liability scale and vice versa, so weight qualitative (sign-changing) interactions over removable quantitative ones, and check robustness across scales before claiming biology.

Want me to sketch the REGENIE interaction command plus the per-stratum-GWAS → cross-trait LDSC pipeline as a runnable GWF/pixi workflow for GenomeDK?

## GWF Workflow

I'll build the workflow files, then walk you through running it. Let me check the helper scripts compile cleanly before handing them over.Now the syntax check:Let me functionally test the data-handling scripts with toy inputs, since those carry the logic:All scripts run correctly: strata coding handles case, the `PC×Hap` products zero out for R as expected, effect allele maps to `ALLELE1` for munge, and the interaction filter picks the right rows.Keep the `scripts/` files in a `scripts/` subdir next to `workflow.py` (the workflow references `{ROOT}/scripts/...`).

**Run order.** Edit the CONFIG block, then `gwf -f workflow.py run`. The graph is `qc` + `strata` + `int_covars` → `step1_full` → `interaction` → `top_int` for Arm A, and per stratum `gwas_{I,R}` → `munge_{I,R}` → `h2_{I,R}` plus `rg` for Arm B.

A few decisions I made that are worth a second look before you launch:

The interaction covar file carries `Hap` plus `PC1_x_Hap … PC20_x_Hap`, and step 2 uses `--interaction Hap`. That gives you the Keller E×covariate guard. What REGENIE won't add is the per-SNP G×PC term, so treat any genome-wide-significant interaction hit as provisional until you re-test it inside an ancestry-homogeneous subset or with local ancestry — that's where a real I-vs-R hit and an ancestry artifact diverge.

I defaulted `STRATUM_SPECIFIC_STEP1 = False`, so both per-stratum GWASs reuse the full-sample step-1 predictors. It's cheaper and biases r_g slightly toward 1 (conservative for your hypothesis). If compute allows, flip it to `True` for the cleaner estimate. Either way, with the cohort halved, the LDSC r_g and h2 standard errors are the binding constraint — below roughly 5k per stratum the rg won't be informative, and that's the number to sanity-check first.

Two things I left as deliberate stubs: LDSC runs on autosomes only, so the X chromosome (your most mechanistically interesting target, given the X-linked regulators) needs its own REGENIE step-2 pass with male hemizygous coding. And the binary `Hap` is a coarse stand-in for the heterochromatin-sink variable — if you derive a continuous Yq12 read-depth measure, you can swap it in wherever `Hap` appears and get a more direct test of the dosage model, though `--interaction` on a continuous covariate changes the per-SNP interpretation to a slope.

Want me to add the X-chromosome step-2 target with proper male coding, or wire in the two-step interaction screen (marginal-variance filter → interaction test on survivors) to recover power on the genome-wide pass?




---- 


> In @notes.md (accompanying the initial workflow.py) you wrote "What REGENIE won't add is the per-SNP G×PC term, so treat any genome-wide-significant interaction hit as provisional until you re-test it inside an ancestry-homogeneous subset or with local ancestry — that's where a real I-vs-R hit and an ancestry artifact diverge.". I am not sure I understand why that is not included. Is the SNPxPC interaction not compeltely standard?
Read notes.md (unchanged)

Good question — and it gets at a subtle but important distinction. Let me first confirm exactly what REGENIE's --interaction does, so I'm not hand-waving.

Confirmed: --interaction takes one variable name, and REGENIE fits Y ~ g + Hap + g×Hap + covariates where the covariates (PCs, PC×Hap, age, batch) enter only as additive main effects. It never crosses the tested SNP g with the PCs. So the SNP×PC terms simply aren't part of what that routine estimates.

Now to your actual question — and you're right to be suspicious, because there are two different "PC interactions" and only one of them is the easy, standard one.

**The asymmetry: per-person-constant vs per-SNP**

- PC×Hap (what we did include): PC_i × Hap depends only on the person, not on which SNP is being tested. So you can precompute all 20 columns once and drop them into the covariate file as ordinary main-effect covariates. Cheap, and REGENIE happily carries them.
- SNP×PC (what's missing): g × PC_i changes for every variant in the scan. It can't be precomputed into a covariate file — it has to be fit inside the per-variant model, alongside g and g×Hap. That's a structurally different thing, and REGENIE's single---interaction machinery (with its Firth/SPA test on one interaction parameter) doesn't do it.

So it's not that SNP×PC is exotic — it's that one is a static covariate and the other is a per-SNP model term, and --interaction only handles a single per-SNP interaction.

**Why both terms matter here (the Keller point)**

Keller (2014) says: to test G×E cleanly when a covariate C confounds it, you need both C×E and C×G in the model, not just C as a main effect. With C = PCs (ancestry) and E = Hap:

- PC×Hap (C×E) closes the channel where ancestry modifies the Hap main effect.
- SNP×PC (C×G) closes the channel where ancestry modifies the SNP effect.

Main-effect PCs only control ancestry in the mean phenotype; PC×Hap controls ancestry in the Hap effect; neither touches ancestry in the SNP effect. That last one is the open door, and here it's wide open because Hap (I vs R) is essentially a proxy for the ancestry gradient.

Concretely the artifact looks like this: suppose a SNP's true effect genuinely differs across the Steppe↔WHG gradient — easily produced by differential LD (your genotyped SNP tags the causal variant differently in R-ancestry vs I-ancestry backgrounds), no real epistasis needed. That is a SNP×ancestry interaction. Because the model can only attribute effect-modification to Hap, not to continuous ancestry, it loads that SNP×ancestry signal onto SNP×Hap → a significant, entirely spurious "interaction." A SNP×PC term would let the model assign the modification to ancestry directly, leaving SNP×Hap clean.

**Why it's still not "just add SNP×PC and move on"**

Even if you did fit g×PC1..20 per SNP, global PCs capture global ancestry, while the confounder that bites a specific locus is local ancestry at that locus — which global PCs only smear. So SNP×PC is a partial fix at best. That's why the notes point to the two things that actually remove both global and local ancestry confounding:

- re-test hits in an ancestry-homogeneous subset (PCs ~constant → SNP×PC can't vary → if SNP×Hap survives, it isn't ancestry), or
- local ancestry at the locus as the adjusting/interacting term.

A real I-vs-R epistasis hit persists under both; an ancestry artifact collapses. That divergence is the whole test.

So: SNP×PC as a concept is standard and exactly the right guard in principle; it's just (a) not something REGENIE's --interaction fits, and (b) an incomplete fix versus the homogeneous-subset / local-ancestry check, which is strictly better. That's why I left it as a post-hoc step rather than baking a half-measure into the genome-wide pass.

If useful, I can add a small post-hoc helper that takes the genome-wide-significant hits and re-fits each one with an augmented model (g + Hap + g×Hap + PCs + Hap×PCs + g×PCs) and/or restricted to an ancestry-homogeneous subsample, reporting whether the SNP×Hap term survives. Want me to wire that in?


❯ What do you consider a ancestry-homogeneous subset?

Two things can be meant by "homogeneous," and the distinction is the whole point:

- (a) Everyone ancestrally similar — low total ancestry variance (a narrow slab of PC space).
- (b) Hap decoupled from ancestry — within the subset, the I and R groups have the same ancestry distribution, i.e. Hap ⫫ ancestry.

(b) is what you actually need, and it's weaker (so less costly in N) than (a). The confound is Hap being a proxy for ancestry; once Hap no longer tracks ancestry inside the subset, a surviving SNP×Hap signal can't be an ancestry artifact — regardless of how much ancestry variation remains. (a) achieves (b) trivially by killing the variance, but it throws away far more samples than necessary, and you're already paying a ~4× power penalty on interactions.

**Operational definition**

Ancestry-homogeneous with respect to the I/R contrast = PCs no longer predict Hap. That's directly checkable: fit Hap ~ PC1..PCk (logistic) in the subset and require the PCs to be jointly near-uninformative — AUC ≈ 0.5, pseudo-R² ≈ 0, joint test non-significant. That metric is the definition, not a proxy for it.

**How to get there (largest subset that satisfies it)**

1. Common support first. Drop the PC-space tails where only R (Steppe-rich) or only I (WHG-rich) exist — there's no one to compare against there anyway. Keep the overlap region.
2. Match I to R on the top PCs within that overlap: nearest-neighbour/caliper matching on PC1..PCk, or propensity-score matching where the propensity is P(Hap=I | PCs). This balances the ancestry distribution between groups by construction and keeps the most samples.
3. Verify with the Hap ~ PCs check above. Tighten the caliper until PCs stop separating I from R.

A blunt alternative is just a narrow PC band around the overlap centroid (approach (a)) — simplest, but sheds the most N.

**The caveat that limits all of this**

Global PCs capture global ancestry. The confounder that bites a specific hit is local ancestry at that locus, which PCs only smear. So a PC-matched subset reduces but doesn't fully guarantee removal of local-ancestry confounding at one SNP. The strictly cleaner re-test conditions on (or restricts to common) local ancestry at the hit locus — RFMix/FLARE-style local-ancestry calls. PC-matching is the cheap, genome-wide-applicable approximation; local ancestry is the gold standard for a handful of top hits.

So in practice: for the genome-wide pass, trust nothing; for the surviving hits, re-fit each one in a PC-matched I/R subset (verified Hap ⫫ PCs), and for the few that still hold up, confirm with local ancestry at the locus. A real I-vs-R epistasis hit survives all three; an artifact dies at the first or second.

If you want, I can add a helper that builds the matched subset (PC caliper or propensity matching), prints the Hap ~ PCs AUC before/after as the homogeneity check, and emits a --keep list to re-run the top hits through REGENIE --interaction. Want that?


❯ I guess this whole exercise because I don't have that. I only have the data I have described. Would ancestry-matched permutation be a way forward?

Short answer: yes — ancestry-matched permutation is a sound and well-judged move with exactly the data you have, and it's strictly better than either trusting the asymptotic p-values or doing naive label permutation. But it buys you protection only up to the resolution of your global PCs, and it's worth being precise about why it works and where it stops.

**Why it works (and why naive permutation doesn't)**

The confound is that Hap is a proxy for ancestry, so SNP×Hap absorbs SNP×ancestry. The two permutation schemes differ in what they preserve:

- Free permutation of Hap destroys both the Hap↔ancestry correlation and any true epistasis. Beating that null tells you nothing about ancestry — it's the wrong null.
- Permute Hap only within ancestry-matched sets (k-NN neighbourhoods or fine PC strata). This keeps the Hap↔ancestry relationship intact (you only swap labels between people of the same ancestry) while breaking the association between Hap and genotype beyond what ancestry implies. Within a stratum where ancestry is ~constant, the SNP genotype distribution is the same for I and R — so under the permutation, the SNP×Hap statistic has no confounding and no epistasis. The observed statistic exceeding that null is evidence of real interaction.

So it's a non-parametric implementation of the matched-subset idea, but it uses the whole sample (every stratum contributes), which is much better for power than carving out one homogeneous subset.

**The design, concretely**

1. Build ancestry-matched sets: k-nearest-neighbours on PC1..PCk, or fine PC clusters with both I and R present.
2. Permute Hap within each set; hold genotype, phenotype, and all covariates fixed.
3. Recompute the SNP×Hap statistic from the full model (g + Hap + g×Hap + PCs + Hap×PCs) on the permuted labels.
4. Empirical p = fraction of permutations with statistic ≥ observed.

**Two gotchas worth respecting:**

- Keep both main effects in the model. Permuting Hap also breaks its main effect on autism; if you test the bare interaction coefficient from a model that already includes Hap and g main effects, the statistic isolates the interaction and the calibration is approximately valid. (If you want it exact in the presence of strong main effects, permute model residuals instead — but for a screen this is overkill.)
- Stratum fineness = your control. Coarse strata leave residual within-stratum ancestry variation for Hap to proxy. Tighten until Hap ~ PCs within strata is uninformative.

A nice cheap by-product: run a modest number of genome-wide permutations and compare the observed interaction λ_GC to the permutation-null λ. That's a global diagnostic for how much of your scan-wide inflation is just Hap↔ancestry structure — far more honest than comparing λ to 1. Then spend the expensive per-SNP empirical permutation only on the handful of top hits.

**The hard limit — say it out loud**

Permutation cannot manufacture information you don't have. It controls confounding to the resolution of global PCs, exactly like the matched subset. If at a specific locus the local ancestry differs between I and R beyond what global PCs capture (entirely plausible near loci under differential selection), then within every global-PC stratum I and R still differ in genotype there, and the ancestry-matched permutation will not null it out — the artifact sails through and looks real. So a survivor means "not explained by global ancestry structure," not "proven epistasis."

That residual is irreducible with autosomal genotypes + global PCs alone. The things that actually close it need data you said you don't have: local-ancestry calls at the hit locus, or replication in an independent cohort. So the honest pipeline is: ancestry-matched permutation as the genome-wide-credible filter, and for survivors, flag them explicitly as needing local ancestry or replication before any biological claim.

One genuinely orthogonal lever you do have: your mechanistic prior. A survivor that lands in a pre-specified set (imprinted clusters, X-linked regulators, your meiosis genes) is more credible than a random genome-wide survivor, because the prior is independent of the ancestry-artifact process. Restricting the permutation test to that candidate set also slashes the multiple-testing burden — which, given the 4× interaction power penalty, may matter more than anything else here.

Want me to implement the ancestry-matched permutation? Concretely: a helper that builds PC-kNN matched sets, runs within-set Hap permutations, emits both the genome-wide λ-vs-null diagnostic and per-hit empirical p-values for the top interactions — wired in as a perm_interaction target after top_int.


---

The perm_interaction target — ancestry-matched permutation re-test

This is the G-side guard REGENIE's --interaction can't provide, wired in as a post-hoc target after top_int. scripts/ancestry_matched_perm.py:

1. Picks the top interaction hits (+ a random SNP panel) from the REGENIE output, reads their genotypes straight from the .bed (a small built-in SNP-major reader — no plink round-trip).
2. Stratifies on the Hap propensity score P(Hap|PCs) and permutes Hap within strata. A real I-vs-R effect lives within strata → broken by permutation → survives; an ancestry artifact lives between strata → preserved → collapses.
3. Auto-selects strata granularity by your tightening rule: make strata fine enough that Hap is unpredictable from PCs within a stratum (AUC→0.5), then take the loosest such for power. The selection curve is written to perm_strata_selection.tsv.
4. Emits per-hit anc_matched_emp_p, plus a genome-wide λ-vs-permutation-null diagnostic.

The principled tightening — and a real discovery

Your instinct to ask for it was well-placed: my first cut (kmeans on all 10 PCs) never decoupled (within-stratum AUC stuck at 0.63–0.75), because Euclidean matching wastes resolution on the 9 noise PCs. The fix is textbook and is exactly the "tighten until Hap ~ PCs is uninformative" criterion made operational: stratify on the 1-D Hap-propensity score (the only direction that carries the confounding). That decouples cleanly (AUC ~0.49) and lets the loosest granularity pass, maximizing power. The AUC-vs-granularity curve is the selection rule, reported for audit.

Proven correct, not just runnable

The synthetic data now plants a confound-only SNP (effect varies with ancestry, no true Hap interaction) alongside the real interaction. In the test they have near-identical naive interaction statistics (23.4 vs 21.5) — indistinguishable to the scan — yet the permutation cleanly separates them:

- real rs3904: anc_matched_emp_p = 0.003 (survives)
- confound rs3779: anc_matched_emp_p = 0.375 (collapses)

That separation is asserted in the test suite, which is fully green end-to-end.

The honest limit is documented everywhere (README caveat 1, CLAUDE.md, the script header): this controls confounding only to global-PC resolution — local ancestry at a locus survives it, so survivors remain provisional pending local-ancestry re-test or replication.

All PERM_* knobs are env-overridable (defaults: top 500, 1000 perms/hit); the target is in the default graph. Nothing committed, per your git defaults.