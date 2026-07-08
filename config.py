"""Central configuration for the Bayesian oracle (2A/2B information ceiling).

Every path and numeric knob used by the oracle scripts lives here, so the rest of the
code has no hard-coded constants. All values are **environment-overridable**: set the
matching `ORACLE_*` variable (or the legacy `MAX_2B_N`) to point the oracle at your own
data / Python-2 env / scoring harness without editing any source. Defaults reproduce the
paper's run on this machine.

    from config import PATHS, PARAMS, sample_list

Paper reference: `phasing_gate0_figures/article.tex`, Prop. `prop:oraclebound`,
Cor. `cor:decodegap`, Eq. `eq:post2b` / `eq:vaf`.
"""
import os

# --------------------------------------------------------------------------------------
# Layout
# --------------------------------------------------------------------------------------
# ROOT defaults to the parent of this file's directory (i.e. DECODE_2_3), matching the
# original `dirname(dirname(abspath(__file__)))`. Override ORACLE_ROOT for a relocated
# checkout whose data live elsewhere.
_HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.environ.get("ORACLE_ROOT", os.path.dirname(_HERE))
PARENT = os.environ.get("ORACLE_PARENT", os.path.dirname(ROOT))


def _p(env, *parts, base=ROOT):
    """Env-overridable path: ORACLE_<env> wins, else base/<parts...>."""
    return os.environ.get(f"ORACLE_{env}", os.path.join(base, *parts))


PATHS = {
    # Python-2 interpreter of the SMC-Het scoring env (least portable; override per host).
    "python2": os.environ.get(
        "ORACLE_PYTHON2",
        "/users/eleves-b/2023/louis.darrigol/.conda/envs/smchet-py2/bin/python",
    ),
    # The SMC-Het scoring harness (same one DECODE uses), imported with PARENT on PYTHONPATH.
    "smcscoring": _p("SMCSCORING", "DECODE2_reviewed", "smc_het_eval", "SMCScoring.py",
                     base=PARENT),
    "parent": PARENT,                                   # goes on PYTHONPATH for the harness
    # Inputs.
    "truth": _p("TRUTH", "scoring", "truth"),           # per-tumour truth.{1A,1C,2A,2B,vcf}
    "dream": _p("DREAM", "DREAM_data"),                 # *_mutation_table_with_multiplicity.csv
    # Outputs / comparison.
    "oracle_scores": _p("OUT", "scoring", "oracle_scores"),   # oracle 2A/2B written here
    "decode_scores": _p("SCORES", "scoring", "scores"),       # DECODE 2A/2B for the comparison
}

# --------------------------------------------------------------------------------------
# Numeric knobs for the oracle likelihood (oracle_score_one.py)
# --------------------------------------------------------------------------------------
PARAMS = {
    # Largest tumour (by #SNVs in the scoring VCF) for which the N*N 2B CCM is built.
    "max_2b_n": int(os.environ.get("MAX_2B_N", os.environ.get("ORACLE_MAX_2B_N", "20000"))),
    # Drop truth clusters with CCF at or below this (numerical / empty clusters).
    "min_ccf": float(os.environ.get("ORACLE_MIN_CCF", "0.01")),
    # Clip the per-locus expected VAF xi into (eps_xi, 1-eps_xi) before the binomial log-lik.
    "eps_xi": float(os.environ.get("ORACLE_EPS_XI", "1e-9")),
    # Additive floor on the mixture prior inside the log (avoid log 0).
    "eps_prior": float(os.environ.get("ORACLE_EPS_PRIOR", "1e-12")),
}


# --------------------------------------------------------------------------------------
# Shared helpers (used by both the SLURM worker and the local runner)
# --------------------------------------------------------------------------------------
def vcf_n(path):
    """Number of records (non-comment lines) in a VCF."""
    return sum(1 for line in open(path) if not line.startswith("#"))


def sample_list(max_n=None):
    """Tumour ids under PATHS['truth'] whose scoring VCF has <= max_n records.

    max_n defaults to PARAMS['max_2b_n']. Ordering is deterministic (sorted), so it
    matches the SLURM array index used by oracle.sbatch.
    """
    if max_n is None:
        max_n = PARAMS["max_2b_n"]
    out = []
    for d in sorted(os.listdir(PATHS["truth"])):
        v = os.path.join(PATHS["truth"], d, f"{d}.truth.scoring_vcf.vcf")
        if os.path.isfile(v) and vcf_n(v) <= max_n:
            out.append(d)
    return out


def parse_final_score(s):
    """Parse an SMC-Het score string ``[[AUPR, AJSD], final]`` -> the final combined score.

    Extracts the last comma-separated field (the combined score). Raises ValueError on a
    malformed string rather than silently returning NaN, so a corrupt score file surfaces
    loudly instead of poisoning the aggregate / the paper table downstream.
    """
    try:
        return float(s.rsplit(",", 1)[-1].strip(" ]\n"))
    except (ValueError, IndexError) as e:
        raise ValueError(f"malformed SMC-Het score string: {s!r}") from e
