"""
lag_check.py -- exploratory: is the phase lag larger under CSD than under the average
reference in ds003702 gamma_high? (Yoshinaga et al., 2020: ImC outperforms wPLI at
lags of about pi/4 or more, wPLI at small lags.)

Reads battery_phasediff.csv from the P2 deriv_dir. For every subject x condition x pair
cell, lag_deg = distance of the mean phase difference from the nearest of 0 and 180 deg
(0 = pure zero-lag, 90 = maximal lag). Cells with resultant_len < --rmin are dropped.

Prints, per decomposition x reference: median lag_deg, fraction of cells with lag >= 45 deg,
median resultant length; then a per-subject paired comparison CSD vs average
(ROI-median lag per subject, Wilcoxon signed-rank).

Usage:  python lag_check.py --p2 "<P2 deriv_dir>" [--band gamma_high] [--rmin 0.1]
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--p2", type=Path, required=True)
ap.add_argument("--band", default="gamma_high")
ap.add_argument("--rmin", type=float, default=0.1)
a = ap.parse_args()

ph = pd.read_csv(a.p2 / "battery_phasediff.csv")
ph = ph[ph.band == a.band].copy()
ang = np.abs(((ph.mean_angle_deg + 180.0) % 360.0) - 180.0)          # 0..180
ph["lag_deg"] = np.minimum(ang, 180.0 - ang)                          # 0..90
n_all = len(ph)
ph = ph[ph.resultant_len >= a.rmin]
print(f"band {a.band}: kept {len(ph)} of {n_all} cells with R >= {a.rmin}\n")

summ = (ph.groupby(["decomposition", "reference"])
          .agg(cells=("lag_deg", "size"),
               median_lag=("lag_deg", "median"),
               frac_lag_ge45=("lag_deg", lambda x: (x >= 45).mean()),
               median_R=("resultant_len", "median"))
          .round(3))
print(summ.to_string(), "\n")

for dec, g in ph.groupby("decomposition"):
    s = g.groupby(["subject", "reference"]).lag_deg.median().unstack("reference").dropna()
    refs = list(s.columns)
    csd = [r for r in refs if r.lower().startswith("csd")][0]
    avg = [r for r in refs if r.lower().startswith("av")][0]
    d = s[csd] - s[avg]
    w = stats.wilcoxon(s[csd], s[avg])
    print(f"{dec}: per-subject median lag  avg {s[avg].median():.1f} deg, csd {s[csd].median():.1f} deg, "
          f"diff median {d.median():+.1f} deg, subjects with larger lag under CSD {int((d > 0).sum())}/{len(d)}, "
          f"Wilcoxon p = {w.pvalue:.4g}")
