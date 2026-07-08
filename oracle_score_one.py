#!/usr/bin/env python3
"""
SLURM per-tumour Bayes-oracle 2A/2B scorer (the achievable information ceiling).
Driven by SLURM_ARRAY_TASK_ID over the tumours with vcf_n <= MAX_2B_N. Computes the oracle soft
assignment over the TRUE clusters (knows truth.1C CCFs; Bayes-optimal per read count), writes the N*N
CCM to node-local /tmp (off the 30 GB home quota), scores 2A+2B with the SMC harness, and writes only
the small score to scoring/oracle_scores/<sid>.txt.
"""
import os, sys, gzip, subprocess, csv, shutil
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import PATHS, PARAMS, sample_list, parse_final_score

PYTHON2 = PATHS["python2"]
SMCSCORING = PATHS["smcscoring"]
PARENT = PATHS["parent"]
TRUTH = PATHS["truth"]
DREAM = PATHS["dream"]
OUT = PATHS["oracle_scores"]


def parse_vcf(p):
    rows = []
    for line in open(p):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        ad = f[10].split(":")[1].split(",")
        rows.append((f[0], int(f[1]), int(ad[0]), int(ad[1])))
    return rows


def parse_1c(p):
    return [(int(f[1]), float(f[2])) for f in (l.split() for l in open(p)) if len(f) >= 3]


def load_cnm(sid):
    fp = os.path.join(DREAM, f"{sid}_mutation_table_with_multiplicity.csv")
    d = {}
    if os.path.isfile(fp):
        for r in csv.DictReader(open(fp)):
            d[(r["Chromosome"], int(r["Position"]))] = (
                int(r["major_cn"]) + int(r["minor_cn"]), float(r["multiplicity"]))
    return d


def oracle_R(sid):
    """Bayes-oracle responsibilities R[i,k] = P(mutation i belongs to true cluster k | read count).

    This is the posterior-mean predictor of Eq. (eq:post2b), evaluated at the TRUE cluster
    prevalences/depths -- i.e. the numerical realisation of the closed-form 2B ceiling S_2B*
    (Prop. prop:oraclebound). It suffers only read noise + the equal-CCF degeneracy.
    """
    d = os.path.join(TRUTH, sid)
    vcf = parse_vcf(os.path.join(d, f"{sid}.truth.scoring_vcf.vcf"))
    # True clusters: CCF phi_k (truth.1C, drop empty ones), purity rho (truth.1A), copy number.
    cl = [c for c in parse_1c(os.path.join(d, f"{sid}.truth.1C.txt")) if c[1] > PARAMS["min_ccf"]]
    rho = float(open(os.path.join(d, f"{sid}.truth.1A.txt")).read().split()[0])
    cnm = load_cnm(sid)
    # ccf = phi_k (true prevalences); pi = mixture prior over clusters ∝ #mutations, sum to 1.
    ccf = np.array([c[1] for c in cl]); pi = np.array([c[0] for c in cl], float); pi /= pi.sum()
    R = np.zeros((len(vcf), len(cl)))
    for i, (ch, pos, ref, alt) in enumerate(vcf):
        b, dd = alt, ref + alt                       # b = variant reads, dd = depth d_i = ref+alt
        C, m = cnm.get((ch, pos), (2, 1.0))          # locus copy number C_n and multiplicity m_n
        eps = PARAMS["eps_xi"]
        # Expected VAF under each cluster, Eq. (eq:vaf): xi_k = phi_k * m / (rho*C + 2(1-rho)).
        xi = np.clip(ccf * m / (rho * C + 2 * (1 - rho)), eps, 1 - eps)
        # Log-posterior over clusters: log prior + binomial log-likelihood log Binom(b | d, xi_k),
        # i.e. the per-locus factor of the read-count model (Eq. eq:lik). Binomial coeff drops out
        # (constant in k, killed by the softmax normalisation below).
        lp = np.log(pi + PARAMS["eps_prior"]) + b * np.log(xi) + (dd - b) * np.log1p(-xi)
        # Softmax -> normalised posterior responsibilities (subtract max for numerical stability).
        lp -= lp.max(); r = np.exp(lp); R[i] = r / r.sum()
    return R


def main():
    idx = int(os.environ.get("SLURM_ARRAY_TASK_ID", sys.argv[1] if len(sys.argv) > 1 else "1"))
    sids = sample_list()
    if idx < 1 or idx > len(sids):
        print(f"idx {idx} out of range 1..{len(sids)}"); return
    sid = sids[idx - 1]
    print(f"### {sid} (idx {idx}/{len(sids)}) node={os.uname()[1]}", flush=True)
    R = oracle_R(sid)
    d = os.path.join(TRUTH, sid)
    vcf = os.path.join(d, f"{sid}.truth.scoring_vcf.vcf")
    env = dict(os.environ); env["PYTHONPATH"] = PARENT + os.pathsep + env.get("PYTHONPATH", "")
    tmp = os.path.join(os.environ.get("TMPDIR", "/tmp"), f"oracle_{sid}_{os.getpid()}")
    os.makedirs(tmp, exist_ok=True)
    res = {}
    # try/finally guarantees the (large, off-quota) node-local scratch is removed even if
    # scoring raises; the exception itself is NOT swallowed: it propagates so the SLURM
    # array task exits non-zero rather than writing a poisoned score.
    try:
        # 2A = hard assignment: MAP cluster per mutation (argmax posterior, 1-indexed).
        p2a = os.path.join(tmp, f"{sid}_estimate.2A.txt")
        open(p2a, "w").write("\n".join(str(int(x)) for x in R.argmax(1) + 1) + "\n")
        # 2B = soft co-clustering matrix: CCM[i,j] = P(z_i = z_j | data) = sum_k R[i,k] R[j,k]
        # (posterior prob. that i,j share a cluster). Diagonal forced to 1 (i co-clusters with itself).
        p2b = os.path.join(tmp, f"{sid}_estimate.2B.txt.gz")
        ccm = R @ R.T; np.fill_diagonal(ccm, 1.0)
        with gzip.open(p2b, "wt") as fh:
            np.savetxt(fh, ccm, fmt="%.6g", delimiter="\t")
        del ccm
        for ch, pred, tf in [("2A", p2a, f"{sid}.truth.2A.txt"), ("2B", p2b, f"{sid}.truth.2B.gz")]:
            o = os.path.join(tmp, f"{sid}.{ch}.txt")
            cmd = [PYTHON2, SMCSCORING, "-c", ch, "--predfiles", pred,
                   "--truthfiles", os.path.join(d, tf), "--vcf", vcf, "-o", o]
            with open(o + ".log", "w") as log:
                rc = subprocess.call(cmd, cwd=PARENT, env=env, stdout=log, stderr=subprocess.STDOUT)
            if rc != 0 or not os.path.exists(o):
                raise RuntimeError(
                    f"{sid}: SMCScoring {ch} failed (rc={rc}, output={'present' if os.path.exists(o) else 'missing'}); "
                    f"see {o + '.log'}")
            raw = open(o).read().strip()
            parse_final_score(raw)          # validate the format now (fail-fast); store the raw string
            res[ch] = raw
    finally:
        shutil.rmtree(tmp, ignore_errors=True)
    os.makedirs(OUT, exist_ok=True)
    with open(os.path.join(OUT, f"{sid}.txt"), "w") as f:
        f.write(f"2A\t{res['2A']}\n2B\t{res['2B']}\n")
    print(f"### {sid} oracle 2A={res['2A']} 2B={res['2B']}", flush=True)


if __name__ == "__main__":
    main()
