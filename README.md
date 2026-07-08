# Bayesian oracle (2A/2B information ceiling)

Isolated scripts for the **Bayes oracle** of the paper
(`phasing_gate0_figures/article.tex`, §Discussion "DECODE already sits at the
information ceiling", and Methods Prop. `prop:oraclebound` / Cor. `cor:decodegap`).

The oracle is the **posterior-mean predictor** of Eq. `eq:post2b`, evaluated at the
**true** cluster prevalences and depths. It is *not* an extra model — it is the
numerical realisation of the closed-form 2B ceiling `S_2B*`, i.e. the best score any
soft reconstruction can attain at fixed depth/data. It suffers only read noise plus
the equal-CCF degeneracy (Prop. `prop:nonid`–`prop:depth`).

## What each mutation gets

For mutation `i` with read count `b_i` at depth `d_i`, over the true clusters `k` with
CCF `c_k` (from `truth.1C`), purity `rho` (`truth.1A`), and per-locus copy number/
multiplicity `(C, m)` (`DREAM_data/*_mutation_table_with_multiplicity.csv`):

    xi_k = c_k * m / (rho*C + 2*(1-rho))          # Eq. (eq:vaf), diploid -> c_k/2
    R[i,k] ∝ pi_k * Binom(b_i | d_i, xi_k)        # posterior responsibility, Eq. (eq:post2b)

The soft co-clustering matrix is `CCM = R @ R.T` (diagonal set to 1); the hard 2A
assignment is `argmax_k R[i,k]`. Both are scored with the **same** SMC-Het harness
DECODE uses (`DECODE2_reviewed/smc_het_eval/SMCScoring.py`, py2 env).

## Files

| file | role |
|------|------|
| `config.py` | **Single source of config.** All paths and numeric knobs, every one env-overridable (see below). Exposes `PATHS`, `PARAMS`, and the shared `sample_list()` / `vcf_n()` helpers. The other scripts hold no hard-coded constants. |
| `oracle_score_one.py` | SLURM per-tumour worker. `SLURM_ARRAY_TASK_ID` -> one tumour (N<=`max_2b_n`). Builds `R`, the N×N CCM on node-local `$TMPDIR`, scores 2A+2B, writes the small score to `oracle_scores/<sid>.txt`. |
| `oracle.sbatch` | SLURM array launcher (`--array=1-45`). Self-locating (works from any checkout); runs `oracle_score_one.py`. |
| `oracle_aggregate.py` | Reads the oracle scores vs DECODE scores, prints per-tumour deltas sorted by 2B gap and the headline means — the paper's **oracle 0.66 vs DECODE 0.62** over the 45 tumours with N<=20000. |

`.oracle_jid` holds the last submitted SLURM job id.

## Configuration (`config.py`)

Nothing needs editing to reproduce the paper on this host. To point the oracle at a
different checkout / data / scoring env, set env vars (defaults in parentheses):

| variable | meaning |
|----------|---------|
| `ORACLE_ROOT` | repo/data root; everything below is relative to it (parent of `bayesian_oracle/`) |
| `ORACLE_PYTHON2` | Python-2 interpreter of the SMC-Het scoring env |
| `ORACLE_SMCSCORING` | path to `SMCScoring.py` |
| `ORACLE_TRUTH` / `ORACLE_DREAM` | truth files / mutation-multiplicity tables (`<ROOT>/scoring/truth`, `<ROOT>/DREAM_data`) |
| `ORACLE_OUT` / `ORACLE_SCORES` | oracle output / DECODE scores for comparison (`<ROOT>/scoring/{oracle_scores,scores}`) |
| `MAX_2B_N` | largest tumour (by #SNVs) for which the N×N 2B CCM is built (`20000`) |
| `ORACLE_MIN_CCF` / `ORACLE_EPS_XI` / `ORACLE_EPS_PRIOR` | numeric knobs (`0.01` / `1e-9` / `1e-12`) |

## Run

```bash
# from the repo root (parent of bayesian_oracle/)
sbatch bayesian_oracle/oracle.sbatch            # 45-tumour array -> scoring/oracle_scores/
python3 bayesian_oracle/oracle_aggregate.py     # headline means + per-tumour deltas
# quick local check on one tumour, no SLURM (index into sample_list(); writes its score file):
python3 bayesian_oracle/oracle_score_one.py 1
```

## Downstream

`phasing_gate0_figures/generate_resolution_table.py` reads
`scoring/oracle_scores/` to fill the `2B_{o/d}` column of Table `app:restable`.
The scripts locate everything relative to `ROOT = DECODE_2_3` (two levels up from the
script), so this directory can live anywhere one level under `DECODE_2_3` without
breaking paths.

**Note.** The `dec2B` numbers the aggregate compares against come from DECODE's own 2B
scoring (`tmp_tests/score_2b_one.py`, scoring DECODE's `results/*_dpbbvi.npz`), which is
a separate DECODE-side step, not part of the oracle.
