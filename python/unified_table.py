"""
unified_table.py -- contrasts and phase-difference summary for the unified battery.

Reads, from each dataset folder given (the repository's cfg.deriv_dir):
  battery_persubject.csv  subject, condition, band, measure, decomposition, reference, pair, value
  battery_phasediff.csv   subject, condition, band, decomposition, reference, pair,
                          mean_angle_deg, resultant_len

Writes to --out:
  unified_contrasts.csv   one row per dataset x band x measure x decomposition x reference x
                          level (roi | pair).  ROI rows: full util_contrast_stats output on the
                          per-subject ROI-average difference.  Pair rows: parametric stats plus
                          a max-|t| sign-flip p_FWER across the pairs of that cell (same flips
                          for every pair; exact if n <= 22, else 100,000 flips, seed 1).
  phase_diff_summary.csv  per dataset x band x decomposition x reference x pair, over
                          subject x condition cells: n_cells, n and fraction within +/-20 deg
                          of 180, of 0, of either (chance for 'either' = 80/360 = 0.222),
                          median resultant length.
  anchor_checks.txt       reproduction of already-reported values; any FAIL means a label,
                          sign or pairing problem and the table must not be used.

Contrasts (difference vector per subject):
  P1 (ds004117)  load      = high - low
  P2 (ds003702)  validity  = mean(valid cells) - mean(invalid cells), averaged over cue,
                 the convention of the published preprint (anchor PLV ROI t = -3.92).

Usage
  python unified_table.py --p1 <P1 deriv_dir> --p2 <P2 deriv_dir> --out <dir> [--exclude sub-024]
Either dataset may be omitted.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd

import util_contrast_stats as ucs

CELL = ["band", "measure", "decomposition", "reference"]
PH_CELL = ["band", "decomposition", "reference"]
SPEC = {
    "P1": dict(dataset="ds004117", contrast="load (high - low)",
               bands=["theta", "beta", "beta_low", "beta_high"]),
    "P2": dict(dataset="ds003702", contrast="validity (valid - invalid)",
               bands=["gamma_low", "gamma_high"]),
}
MEASURES = ["PLV", "dwPLI", "ImC"]


# ----------------------------------------------------------------------------- loading
def _norm(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    for c in df.columns:
        if df[c].dtype == object:
            df[c] = df[c].astype(str).str.strip()
    return df


def load(folder: Path, exclude: list[str]):
    v = _norm(pd.read_csv(folder / "battery_persubject.csv"))
    ph = _norm(pd.read_csv(folder / "battery_phasediff.csv"))
    if exclude:
        v, ph = v[~v.subject.isin(exclude)], ph[~ph.subject.isin(exclude)]
    key = ["subject", "condition", "pair"] + CELL
    dup = v.duplicated(key).sum()
    if dup:
        raise ValueError(f"{folder}: {dup} duplicated rows in battery_persubject.csv "
                         "(a resumed run appended a subject twice?)")
    return v, ph


def validity_of(label: str) -> str:
    s = label.lower()
    if re.search(r"invalid|incongruent|\binv\b", s):
        return "invalid"
    if re.search(r"valid|congruent", s):
        return "valid"
    raise ValueError(f"cannot read validity from condition label {label!r}")


def diff_wide(v: pd.DataFrame, ds: str, flip: bool) -> pd.DataFrame:
    """Per-subject difference, one row per subject x cell x pair."""
    idx = ["subject", "pair"] + CELL
    if ds == "P1":
        labs = set(v.condition.unique())
        if labs != {"low", "high"}:
            raise ValueError(f"P1 conditions must be low/high, got {sorted(labs)}")
        w = v.pivot_table(index=idx, columns="condition", values="value", aggfunc="first")
        d = w["high"] - w["low"]
    else:
        v = v.assign(val=v.condition.map(validity_of))
        n_per = v.groupby("val").condition.nunique()
        if set(n_per.index) != {"valid", "invalid"} or (n_per != 2).any():
            raise ValueError(f"P2 needs 2 valid and 2 invalid cells, got {n_per.to_dict()}")
        w = v.groupby(idx + ["val"]).value.mean().unstack("val")   # mean over cue
        d = w["valid"] - w["invalid"]   # sign convention of the published preprint
        if flip:
            d = -d
    return d.rename("diff").reset_index().dropna(subset=["diff"])


# ------------------------------------------------------------------------------- stats
def maxstat_fwer(X: np.ndarray) -> np.ndarray:
    """X: subjects x pairs.  Max-|t| sign-flip FWER p per pair."""
    n, k = X.shape
    sx2 = (X * X).sum(0)
    t_obs = np.abs(ucs._t_from_sum(X.sum(0), sx2, n))
    thr = t_obs - 1e-10 * np.maximum(1.0, t_obs)
    count = np.zeros(k)

    def tally(signs):
        tmax = np.nanmax(np.abs(ucs._t_from_sum(signs @ X, sx2, n)), axis=1)
        return (tmax[:, None] >= thr[None, :]).sum(0)

    if n <= ucs.EXACT_MAX_N:
        total, bits = 1 << n, np.arange(n, dtype=np.int64)
        for s in range(0, total, 1 << 17):
            kk = np.arange(s, min(s + (1 << 17), total), dtype=np.int64)
            count += tally(1.0 - 2.0 * ((kk[:, None] >> bits) & 1))
        return count / total
    rng = np.random.default_rng(ucs.SEED)
    done = 0
    while done < ucs.N_FLIPS:
        b = min(10_000, ucs.N_FLIPS - done)
        count += tally(rng.choice((-1.0, 1.0), size=(b, n)))
        done += b
    return (count + 1) / (ucs.N_FLIPS + 1)


def contrasts(d: pd.DataFrame, ds: str) -> pd.DataFrame:
    rows = []
    for cell, g in d.groupby(CELL, sort=True):
        base = dict(zip(CELL, cell), dataset=SPEC[ds]["dataset"], contrast=SPEC[ds]["contrast"])
        roi = g.groupby("subject")["diff"].mean()
        npairs = g.groupby("subject").pair.nunique()
        if npairs.nunique() != 1:
            raise ValueError(f"{ds} {cell}: unequal pair counts across subjects")
        rows.append({**base, "level": "roi", "pair": f"ROI({npairs.iloc[0]})",
                     **ucs.contrast_stats(roi.values)})
        X = g.pivot(index="subject", columns="pair", values="diff").dropna()
        fwer = maxstat_fwer(X.values)
        for j, p in enumerate(X.columns):
            r = ucs.summary_from_t(*_t_n(X.values[:, j]))
            rows.append({**base, "level": "pair", "pair": p, **r,
                         "p_fwer": fwer[j], "perm_method": "maxstat"})
    out = pd.DataFrame(rows)
    lead = ["dataset", "contrast", "level"] + CELL + ["pair", "n", "dz", "dz_ci_lo", "dz_ci_hi",
            "t", "df", "p", "p_perm", "p_fwer", "perm_method", "bf10", "bf01", "p_tost"]
    return out[[c for c in lead if c in out.columns]]


def _t_n(x):
    x = x[np.isfinite(x)]
    return float(x.mean() / x.std(ddof=1) * np.sqrt(x.size)), int(x.size)


def phase_summary(ph: pd.DataFrame, ds: str) -> pd.DataFrame:
    a = ((ph.mean_angle_deg + 180.0) % 360.0) - 180.0
    ph = ph.assign(n180=(a.abs() >= 160.0), n0=(a.abs() <= 20.0))
    g = ph.groupby(PH_CELL + ["pair"])
    out = g.agg(n_cells=("n180", "size"), n_near180=("n180", "sum"), n_near0=("n0", "sum"),
                median_R=("resultant_len", "median")).reset_index()
    out["frac_near180"] = out.n_near180 / out.n_cells
    out["frac_near0"] = out.n_near0 / out.n_cells
    out["frac_zero_lag"] = (out.n_near180 + out.n_near0) / out.n_cells
    out["chance_zero_lag"] = 80.0 / 360.0
    out.insert(0, "dataset", SPEC[ds]["dataset"])
    return out


# ----------------------------------------------------------------------------- anchors
def _pick(tab, ds, band, measure, decomp, ref):
    m = tab[(tab.dataset == SPEC[ds]["dataset"]) & (tab.level == "roi") & (tab.band == band)
            & (tab.measure == measure) & (tab.decomposition.str.lower().str.startswith(decomp))
            & (tab.reference.str.lower().str.startswith(ref))]
    return None if m.empty else m.iloc[0]


def _pair_is(p, a, b):
    toks = set(re.split(r"[^A-Za-z0-9]+", str(p)))
    return toks == {a, b}


def anchors(tab, phs, have) -> list[str]:
    lines = []

    def check(label, row, field, target, tol):
        if row is None:
            lines.append(f"MISSING  {label}")
            return
        got = float(row[field])
        ok = abs(got - target) <= tol
        lines.append(f"{'PASS' if ok else 'FAIL':7s}  {label}: {field} = {got:.4f} "
                     f"(target {target}, tol {tol})")

    if "P2" in have:
        check("P2 PLV Hilbert avg gamma_high ROI", _pick(tab, "P2", "gamma_high", "PLV", "hil", "av"),
              "t", -3.92, 0.01)
        check("P2 ImC Hilbert avg gamma_high ROI", _pick(tab, "P2", "gamma_high", "ImC", "hil", "av"),
              "t", -0.146, 0.01)
        check("P2 dwPLI Hilbert avg gamma_high ROI", _pick(tab, "P2", "gamma_high", "dwPLI", "hil", "av"),
              "bf01", 6.18, 0.01)
        check("P2 PLV Hilbert CSD gamma_high ROI", _pick(tab, "P2", "gamma_high", "PLV", "hil", "csd"),
              "t", -0.63, 0.01)
        f = phs[(phs.dataset == "ds003702") & (phs.band == "gamma_high")
                & phs.decomposition.str.lower().str.startswith("hil")
                & phs.reference.str.lower().str.startswith("av")
                & phs.pair.map(lambda p: _pair_is(p, "F2", "P4"))]
        if f.empty:
            lines.append("MISSING  P2 F2-P4 phase difference")
        else:
            r = f.iloc[0]
            ok = (int(r.n_near180), int(r.n_cells)) == (175, 188)
            lines.append(f"{'PASS' if ok else 'FAIL':7s}  P2 F2-P4 Hilbert avg: {int(r.n_near180)}/"
                         f"{int(r.n_cells)} within 20 deg of 180 (target 175/188)")
    if "P1" in have:
        check("P1 dwPLI Hilbert avg beta ROI", _pick(tab, "P1", "beta", "dwPLI", "hil", "av"),
              "bf01", 4.47, 0.01)
        check("P1 dwPLI wavelet avg beta ROI", _pick(tab, "P1", "beta", "dwPLI", "wav", "av"),
              "bf01", 2.17, 0.01)
    return lines


# -------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--p1", type=Path)
    ap.add_argument("--p2", type=Path)
    ap.add_argument("--out", type=Path, default=Path("."))
    ap.add_argument("--exclude", nargs="*", default=[])
    ap.add_argument("--p2-flip", action="store_true")
    a = ap.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)

    tabs, phs, have = [], [], []
    for ds, folder in (("P1", a.p1), ("P2", a.p2)):
        if folder is None:
            continue
        v, ph = load(folder, a.exclude)
        missing = set(SPEC[ds]["bands"]) - set(v.band) or set(MEASURES) - set(v.measure)
        if missing:
            raise ValueError(f"{ds}: missing levels {missing}")
        tabs.append(contrasts(diff_wide(v, ds, a.p2_flip), ds))
        phs.append(phase_summary(ph, ds))
        have.append(ds)
    if not have:
        ap.error("give --p1 and/or --p2")

    tab, phs = pd.concat(tabs, ignore_index=True), pd.concat(phs, ignore_index=True)
    suffix = "" if not a.exclude else "_excl-" + "-".join(a.exclude)
    tab.to_csv(a.out / f"unified_contrasts{suffix}.csv", index=False, float_format="%.6g")
    phs.to_csv(a.out / f"phase_diff_summary{suffix}.csv", index=False, float_format="%.6g")

    lines = anchors(tab, phs, have) if not a.exclude else ["anchors skipped (exclusion run)"]
    (a.out / f"anchor_checks{suffix}.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    roi = tab[tab.level == "roi"]
    with pd.option_context("display.width", 200, "display.max_rows", 200):
        print(roi[["dataset", "band", "measure", "decomposition", "reference", "n", "dz", "t",
                   "p", "p_perm", "bf10", "bf01", "p_tost"]].round(4).to_string(index=False))


if __name__ == "__main__":
    main()
