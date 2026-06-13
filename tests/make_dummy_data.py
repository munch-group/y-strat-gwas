#!/usr/bin/env python3
"""Generate a small, self-consistent synthetic dataset for the GxHaplogroup
pipeline so the whole workflow can be exercised end-to-end on a laptop.

It writes (into --out-dir):
  geno.vcf            multi-chromosome genotypes (imported to plink1 by the runner)
  phenotypes.txt      FID IID autism            (1/0/NA)
  covariates_base.txt FID IID PC1..PCk age batch
  haplogroup.txt      FID IID Hap               (I or R)

The data is deliberately structured to mimic the real analysis problem:

  * Haplogroup I vs R is *confounded with ancestry*: an ancestry axis (captured
    by PC1) differs in mean between the two groups, and SNP allele frequencies
    drift along that axis. This is exactly the structure the PCxHap covariates
    in Arm A are meant to absorb.

  * One designated SNP (printed as TRUE_INTERACTION_SNP) has a genuine
    SNP x Hap effect on the phenotype (effect present in I, absent in R), so a
    correct interaction scan should rank it near the top.

Sizes are tiny (for speed); LDSC h2/rg estimates will be noisy -- the point of
the test is that every stage runs and produces well-formed output, not that the
estimates are precise.
"""
import argparse
import gzip
import os
import numpy as np

N_PER_GROUP = 1000          # samples per haplogroup (I, R)
N_SNP       = 4000          # total autosomal variants, spread across chr 1..22
N_XSNP      = 300           # chrX variants (non-PAR, male hemizygous)
NPC         = 10            # PCs written to the covariate file
SEED        = 12345
# chrX non-PAR window (hg38): between PAR1 end (2,781,479) and PAR2 start
# (155,701,383), so REGENIE codes males hemizygously.
X_POS_START = 3_000_000


def write_ldsc_reference(ref_dir, G, chrom, bp, snp_ids):
    """Emit per-chromosome LDSC reference files (the format --ref-ld-chr/--w-ld-chr
    consume): {chr}.l2.ldscore.gz, {chr}.l2.M, {chr}.l2.M_5_50.

    LD scores are computed directly here (sum of r^2 within the chromosome) so the
    test never needs LDSC's bed reader -- the real pipeline likewise gets reference
    LD scores as input rather than computing them.
    """
    os.makedirs(ref_dir, exist_ok=True)
    for c in range(1, 23):
        idx = np.where(chrom == c)[0]
        Gc = G[:, idx].astype(float)
        p = Gc.mean(0) / 2.0
        sd = Gc.std(0)
        sd[sd == 0] = 1.0
        Z = (Gc - Gc.mean(0)) / sd
        nobs = Z.shape[0]
        r2 = ((Z.T @ Z) / nobs) ** 2                 # |idx| x |idx| correlation^2
        l2 = r2.sum(1)                               # LD score per SNP (incl. self=1)
        maf = np.minimum(p, 1.0 - p)
        with gzip.open("%s/%d.l2.ldscore.gz" % (ref_dir, c), "wt") as f:
            f.write("CHR\tSNP\tBP\tL2\n")
            for k, j in enumerate(idx):
                f.write("%d\t%s\t%d\t%.4f\n" % (c, snp_ids[j], bp[j], l2[k]))
        with open("%s/%d.l2.M" % (ref_dir, c), "w") as f:
            f.write("%d\n" % len(idx))
        with open("%s/%d.l2.M_5_50" % (ref_dir, c), "w") as f:
            f.write("%d\n" % int((maf > 0.05).sum()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--ref-dir", default=None,
                    help="LDSC reference dir (default: <out-dir>/eur_w_ld_chr)")
    a = ap.parse_args()
    os.makedirs(a.out_dir, exist_ok=True)
    ref_dir = a.ref_dir or os.path.join(a.out_dir, "eur_w_ld_chr")
    rng = np.random.default_rng(SEED)

    n = 2 * N_PER_GROUP
    # Hap: first half I (=1), second half R (=0)
    hap = np.array([1] * N_PER_GROUP + [0] * N_PER_GROUP)

    # ---- ancestry axis confounded with haplogroup -------------------------
    # I individuals sit ~+1 SD, R individuals ~-1 SD on a latent ancestry axis;
    # PC1 is a noisy readout of it. PC2.. are pure noise.
    ancestry = rng.normal(loc=np.where(hap == 1, 0.8, -0.8), scale=1.0)
    pcs = rng.normal(size=(n, NPC))
    pcs[:, 0] = ancestry + rng.normal(scale=0.3, size=n)   # PC1 tracks ancestry

    age   = rng.integers(18, 65, size=n).astype(float)
    batch = rng.integers(0, 2, size=n)                     # 2-level categorical

    # ---- variant annotation ----------------------------------------------
    # distribute SNPs across 22 chromosomes, ascending bp within each
    chrom = (np.arange(N_SNP) % 22) + 1
    snp_ids, bp = [], np.zeros(N_SNP, dtype=int)
    per_chr_pos = {c: 0 for c in range(1, 23)}
    for j in range(N_SNP):
        c = int(chrom[j])
        per_chr_pos[c] += int(rng.integers(5000, 20000))
        bp[j] = per_chr_pos[c]
        snp_ids.append("rs%d" % (j + 1))
    base_maf  = rng.uniform(0.05, 0.45, size=N_SNP)
    # each SNP's freq shifts a little along the ancestry axis (structure)
    anc_load  = rng.normal(scale=0.06, size=N_SNP)

    # ---- genotypes: per-person freq = base_maf + anc_load*ancestry --------
    G = np.empty((n, N_SNP), dtype=np.int8)
    for j in range(N_SNP):
        p = np.clip(base_maf[j] + anc_load[j] * ancestry, 0.01, 0.99)
        G[:, j] = rng.binomial(2, p)

    # ---- chrX genotypes: all-male cohort => hemizygous, dosage in {0,2} ----
    x_ids = ["rsX%d" % (j + 1) for j in range(N_XSNP)]
    x_bp  = X_POS_START + np.cumsum(rng.integers(5000, 20000, size=N_XSNP))
    x_maf = rng.uniform(0.05, 0.45, size=N_XSNP)
    Gx = np.empty((n, N_XSNP), dtype=np.int8)
    for j in range(N_XSNP):
        Gx[:, j] = 2 * rng.binomial(1, x_maf[j], size=n)   # carrier => dosage 2
    x_int_idx = int(rng.integers(N_XSNP))                  # planted X interaction

    # ---- phenotype: logistic liability -----------------------------------
    # a few autosomal main-effect SNPs + ONE autosomal and ONE chrX genuine
    # SNP x Hap interaction (effect only in I) + ONE pure-confound SNP whose
    # effect varies along ANCESTRY (not Hap). Because Hap tags ancestry, the
    # confound SNP shows a spurious SNP x Hap signal in the naive scan but should
    # NOT survive the ancestry-matched permutation -- the real ones should.
    main_idx = rng.choice(N_SNP, size=5, replace=False)
    rest = np.setdiff1d(np.arange(N_SNP), main_idx)
    int_idx  = int(rng.choice(rest))
    conf_idx = int(rng.choice(np.setdiff1d(rest, [int_idx])))
    logit = (-0.2
             + 0.30 * G[:, main_idx].sum(axis=1) / np.sqrt(5)
             + 0.50 * G[:, int_idx] * (hap == 1)       # autosomal effect only in I
             + 0.60 * (Gx[:, x_int_idx] / 2) * (hap == 1)  # chrX effect only in I
             + 0.55 * G[:, conf_idx] * ancestry        # effect varies with ANCESTRY (no Hap)
             + 0.15 * ancestry)                        # ancestry affects risk

    # Two LD blocks for the local cross-stratum r_g (LAVA) test:
    #   chr20 = GENUINE divergence: block SNP effects flip sign by Hap (effect
    #           +delta in I, -delta in R), independent of ancestry -> a real
    #           within-ancestry local r_g ~ -1 that should SURVIVE permutation.
    #   chr21 = CONFOUND divergence: block SNP effects are modulated by ANCESTRY
    #           (no Hap term) -> an apparent local r_g divergence that should
    #           COLLAPSE under the ancestry-matched permutation.
    real_block = np.where(chrom == 20)[0]
    conf_block = np.where(chrom == 21)[0]
    delta = rng.normal(0.0, 0.11, size=len(real_block))
    gamma = rng.normal(0.0, 0.06, size=len(conf_block))
    logit += (2 * hap - 1) * (G[:, real_block] @ delta)   # sign flips by Hap
    logit += ancestry * (G[:, conf_block] @ gamma)        # modulated by ancestry
    logit -= logit.mean()
    prob = 1.0 / (1.0 + np.exp(-logit))
    autism = (rng.uniform(size=n) < prob).astype(int)
    # inject a few missing phenotypes
    miss = rng.choice(n, size=max(1, n // 100), replace=False)
    autism_str = autism.astype(object)
    autism_str[miss] = "NA"

    fid = ["FID%d" % i for i in range(n)]
    iid = ["IID%d" % i for i in range(n)]

    # ---- write VCF --------------------------------------------------------
    code = {0: "0/0", 1: "0/1", 2: "1/1"}
    vcf = os.path.join(a.out_dir, "geno.vcf")
    with open(vcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        for c in range(1, 23):
            f.write("##contig=<ID=%d>\n" % c)
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join("%s_%s" % (fid[i], iid[i]) for i in range(n)) + "\n")
        for j in np.lexsort((bp, chrom)):   # VCF must be sorted by (chrom, pos)
            gts = "\t".join(code[int(g)] for g in G[:, j])
            f.write("%d\t%d\t%s\tA\tG\t.\tPASS\t.\tGT\t%s\n"
                    % (chrom[j], bp[j], snp_ids[j], gts))

    # ---- write chrX VCF (hemizygous males encoded as 0/0 or 1/1 = dosage 0/2)
    # plus a sex file so the X fileset's .fam gets SEX=1 (male) for REGENIE.
    xvcf = os.path.join(a.out_dir, "genoX.vcf")
    with open(xvcf, "w") as f:
        f.write("##fileformat=VCFv4.2\n##contig=<ID=X>\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">\n')
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t"
                + "\t".join("%s_%s" % (fid[i], iid[i]) for i in range(n)) + "\n")
        for j in range(N_XSNP):
            gts = "\t".join(code[int(g)] for g in Gx[:, j])
            f.write("X\t%d\t%s\tA\tG\t.\tPASS\t.\tGT\t%s\n"
                    % (int(x_bp[j]), x_ids[j], gts))
    with open(os.path.join(a.out_dir, "sex.txt"), "w") as f:
        for i in range(n):
            f.write("%s\t%s\t1\n" % (fid[i], iid[i]))   # FID IID SEX(1=male)

    # ---- write pheno / covar / hap ---------------------------------------
    def write_table(path, header, rows):
        with open(path, "w") as f:
            f.write(header + "\n")
            for r in rows:
                f.write("\t".join(str(x) for x in r) + "\n")

    write_table(os.path.join(a.out_dir, "phenotypes.txt"),
                "FID IID autism".replace(" ", "\t"),
                [(fid[i], iid[i], autism_str[i]) for i in range(n)])

    pc_names = ["PC%d" % k for k in range(1, NPC + 1)]
    write_table(os.path.join(a.out_dir, "covariates_base.txt"),
                "\t".join(["FID", "IID"] + pc_names + ["age", "batch"]),
                [tuple([fid[i], iid[i]]
                       + ["%.5f" % pcs[i, k] for k in range(NPC)]
                       + [int(age[i]), int(batch[i])]) for i in range(n)])

    write_table(os.path.join(a.out_dir, "haplogroup.txt"),
                "FID IID Hap".replace(" ", "\t"),
                [(fid[i], iid[i], "I" if hap[i] == 1 else "R") for i in range(n)])

    write_ldsc_reference(ref_dir, G, chrom, bp, snp_ids)

    with open(os.path.join(a.out_dir, "truth.txt"), "w") as f:
        f.write("interaction_snp\t%s\n" % snp_ids[int_idx])
        f.write("maineffect_snps\t%s\n" % ",".join(snp_ids[k] for k in sorted(main_idx)))
        f.write("interaction_snp_X\t%s\n" % x_ids[x_int_idx])
        f.write("confound_snp\t%s\n" % snp_ids[conf_idx])
        f.write("real_block_locus\t20:%d-%d\n"
                % (int(bp[real_block].min()), int(bp[real_block].max())))
        f.write("conf_block_locus\t21:%d-%d\n"
                % (int(bp[conf_block].min()), int(bp[conf_block].max())))

    n_case = int((autism[np.array([x != 'NA' for x in autism_str])]).sum())
    print("wrote synthetic data to %s" % a.out_dir)
    print("  samples            : %d (%d I, %d R)" % (n, N_PER_GROUP, N_PER_GROUP))
    print("  variants           : %d across chr 1..22 + %d on chrX" % (N_SNP, N_XSNP))
    print("  cases (non-missing) : %d" % n_case)
    print("  LDSC reference dir  : %s" % ref_dir)
    print("  TRUE_X_INTERACTION  : %s  (chrX, hemizygous)" % x_ids[x_int_idx])
    print("  CONFOUND_SNP        : %s  (g x ancestry, no g x Hap)" % snp_ids[conf_idx])
    print("  TRUE_INTERACTION_SNP: %s  (chr%d)" % (snp_ids[int_idx], chrom[int_idx]))
    print("  TRUE_MAINEFFECT_SNPs: %s" % ",".join(snp_ids[k] for k in sorted(main_idx)))


if __name__ == "__main__":
    main()
