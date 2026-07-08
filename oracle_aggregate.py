#!/usr/bin/env python3
"""Aggregate the SLURM Bayes-oracle run: oracle 2A/2B vs DECODE 2A/2B over the 45 tumours (N<=20000).
Reports per-tumour deltas (sorted by 2B gap) and the headline means."""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from config import PATHS, parse_final_score

ORA = PATHS["oracle_scores"]
DEC = PATHS["decode_scores"]


def read_oracle(sid):
    """(2A, 2B) final scores from an oracle score file. A malformed file raises
    (via parse_final_score) instead of silently degrading to NaN -- if the worker wrote
    this file it must be well-formed. A channel absent from the file is NaN."""
    d = {}
    for line in open(os.path.join(ORA, f"{sid}.txt")):
        ch, _, val = line.partition("\t")
        d[ch.strip()] = parse_final_score(val)
    return d.get("2A", float("nan")), d.get("2B", float("nan"))


def read_decode(sid, ch):
    """DECODE's final score for (sid, ch); NaN if this tumour was never scored by DECODE
    (a legitimate absence), but a present-but-malformed file raises."""
    p = os.path.join(DEC, sid, f"{sid}.{ch}.txt")
    if not os.path.isfile(p):
        return float("nan")
    return parse_final_score(open(p).read())


rows = []
for f in sorted(os.listdir(ORA)):
    if not f.endswith(".txt"):
        continue
    sid = f[:-4]
    o2a, o2b = read_oracle(sid)
    d2a, d2b = read_decode(sid, "2A"), read_decode(sid, "2B")
    rows.append((sid.replace("-noXY", ""), o2a, d2a, o2b, d2b))

rows.sort(key=lambda r: (r[3] - r[4]) if np.isfinite(r[3] - r[4]) else -1, reverse=True)

print(f"{'tumour':<7}{'ora2A':>7}{'dec2A':>7}{'d2A':>7}   {'ora2B':>7}{'dec2B':>7}{'d2B':>7}")
for name, o2a, d2a, o2b, d2b in rows:
    print(f"{name:<7}{o2a:>7.3f}{d2a:>7.3f}{o2a-d2a:>+7.3f}   {o2b:>7.3f}{d2b:>7.3f}{o2b-d2b:>+7.3f}")

A = np.array([(r[1], r[2], r[3], r[4]) for r in rows], float)
ok2a = np.isfinite(A[:, 0]) & np.isfinite(A[:, 1])
ok2b = np.isfinite(A[:, 2]) & np.isfinite(A[:, 3])
print(f"\nn={len(rows)}  (2A valid {ok2a.sum()}, 2B valid {ok2b.sum()})")
print(f"MEAN 2A   oracle={A[ok2a,0].mean():.3f}  DECODE={A[ok2a,1].mean():.3f}  "
      f"gap={A[ok2a,0].mean()-A[ok2a,1].mean():+.3f}")
print(f"MEAN 2B   oracle={A[ok2b,2].mean():.3f}  DECODE={A[ok2b,3].mean():.3f}  "
      f"gap={A[ok2b,2].mean()-A[ok2b,3].mean():+.3f}")
