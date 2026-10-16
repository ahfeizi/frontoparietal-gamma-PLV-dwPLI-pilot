"""
util_contrast_stats.py -- inference for one per-subject difference vector.

Shared verbatim by
  ahfeizi/frontotemporal-WM-pilot              (ds004117, N = 22)
  ahfeizi/frontoparietal-gamma-PLV-dwPLI-pilot (ds003702, N = 47)

contrast_stats(x), x = per-subject differences, returns
  n, mean, sd, dz, t, df, p      one-sample t test against 0, two-sided
  p_perm, perm_method            sign-flip permutation on |t|: exact enumeration of
                                 all 2^n flips if n <= 22, otherwise 100,000 random
                                 flips (seed 1) with p = (k + 1) / (B + 1)
  bf10, bf01                     JZS Bayes factor, two-sided, Cauchy prior scale
                                 r = sqrt(2)/2 = 0.7071 (Rouder et al., 2009); with
                                 r = 0.707 exactly the N = 22 ceiling is 4.4849 and
                                 would print as 4.48, not the reported 4.49
  p_tost                         two one-sided tests, bounds d_z = +/-0.4
  dz_ci_lo, dz_ci_hi             95% CI for d_z by inverting the noncentral t

Run this file directly to execute the self-test against values already reported
in the two preprints, and a brute-force check of the exact sign-flip enumeration.
"""
from __future__ import annotations

import itertools
import math

import numpy as np
from scipy import integrate, optimize, stats

R_SCALE = math.sqrt(2.0) / 2.0   # 'medium' Cauchy scale, 0.7071 (BayesFactor / JASP default)
TOST_BOUND = 0.4
EXACT_MAX_N = 22
N_FLIPS = 100_000
SEED = 1


def jzs_bf10(t: float, n: int, r: float = R_SCALE) -> float:
    """Two-sided JZS BF10 for a one-sample/paired t (Rouder et al., 2009, eq. 1)."""
    nu = n - 1
    base = 1.0 + t * t / nu

    def integrand(g: float) -> float:
        a = 1.0 + n * g * r * r
        return (a ** -0.5
                * ((1.0 + t * t / (a * nu)) / base) ** (-(nu + 1) / 2.0)
                * (2.0 * math.pi) ** -0.5 * g ** -1.5 * math.exp(-1.0 / (2.0 * g)))

    total = 0.0
    for lo, hi in ((0.0, 1.0), (1.0, 100.0), (100.0, math.inf)):
        total += integrate.quad(integrand, lo, hi, limit=500)[0]
    return total


def tost_p(t: float, n: int, bound: float = TOST_BOUND) -> float:
    """TOST p for H0: |d_z| >= bound (bounds expressed in sample-SD units)."""
    df = n - 1
    shift = bound * math.sqrt(n)
    return float(max(stats.t.sf(t + shift, df), stats.t.cdf(t - shift, df)))


def dz_ci(t: float, n: int, level: float = 0.95) -> tuple[float, float]:
    """CI for d_z from the noncentral t distribution."""
    df = n - 1
    a = (1.0 - level) / 2.0

    def f(nc: float, q: float) -> float:
        return float(stats.nct.cdf(t, df, nc)) - q

    def root(q: float) -> float:
        # cdf(t; nc) decreases in nc; widen a bracket around t in unit steps so the
        # solver never evaluates nct far in the tails, where scipy returns NaN.
        lo, hi = t - 1.0, t + 1.0
        while f(lo, q) < 0.0:
            lo -= 1.0
        while f(hi, q) > 0.0:
            hi += 1.0
        return optimize.brentq(f, lo, hi, args=(q,), xtol=1e-10)

    return root(1.0 - a) / math.sqrt(n), root(a) / math.sqrt(n)


def _t_from_sum(S, sx2: float, n: int):
    m = S / n
    var = np.maximum((sx2 - n * m * m) / (n - 1), 0.0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return m / np.sqrt(var / n)


def signflip_p(x, n_flips: int = N_FLIPS, seed: int = SEED,
               exact_max_n: int = EXACT_MAX_N) -> tuple[float, str]:
    """Two-sided sign-flip permutation p on |t|."""
    x = np.asarray(x, dtype=float)
    n = x.size
    sx2 = float(np.sum(x * x))
    t_obs = abs(float(_t_from_sum(x.sum(), sx2, n)))
    thr = t_obs - 1e-10 * max(1.0, t_obs)

    if n <= exact_max_n:
        total = 1 << n
        chunk = 1 << 18
        bits = np.arange(n, dtype=np.int64)
        count = 0
        for start in range(0, total, chunk):
            k = np.arange(start, min(start + chunk, total), dtype=np.int64)
            signs = 1.0 - 2.0 * ((k[:, None] >> bits) & 1)
            count += int(np.count_nonzero(np.abs(_t_from_sum(signs @ x, sx2, n)) >= thr))
        return count / total, "exact"

    rng = np.random.default_rng(seed)
    count, done, chunk = 0, 0, 10_000
    while done < n_flips:
        b = min(chunk, n_flips - done)
        signs = rng.choice((-1.0, 1.0), size=(b, n))
        count += int(np.count_nonzero(np.abs(_t_from_sum(signs @ x, sx2, n)) >= thr))
        done += b
    return (count + 1) / (n_flips + 1), f"mc{n_flips}"


def summary_from_t(t: float, n: int) -> dict:
    """Everything except the permutation p, from t and n alone."""
    df = n - 1
    bf10 = jzs_bf10(t, n)
    lo, hi = dz_ci(t, n)
    return dict(n=n, dz=t / math.sqrt(n), t=t, df=df,
                p=float(2.0 * stats.t.sf(abs(t), df)),
                bf10=bf10, bf01=1.0 / bf10, p_tost=tost_p(t, n),
                dz_ci_lo=lo, dz_ci_hi=hi)


def contrast_stats(x) -> dict:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    n = x.size
    if n < 3:
        raise ValueError(f"need at least 3 finite values, got {n}")
    mean, sd = float(x.mean()), float(x.std(ddof=1))
    if sd == 0.0:
        raise ValueError("zero variance in difference vector")
    t = mean / sd * math.sqrt(n)
    out = dict(n=n, mean=mean, sd=sd)
    out.update(summary_from_t(t, n))
    out["p_perm"], out["perm_method"] = signflip_p(x)
    return out


# --------------------------------------------------------------------------- self-test
# (label, n, t as printed, field, value as printed).  The input t is itself rounded,
# so a check PASSES when the printed value is reachable from some t that rounds to
# the printed t: the field is evaluated at both ends of t's rounding interval and the
# printed value's own rounding interval must overlap that range.  t = 0 entries are
# exact (sample-size ceilings) and test the prior scale directly.
_TARGETS = [
    ("P2 PLV ROI gamma_high", 47, "-3.92", "p", "0.0003"),
    ("P2 PLV ROI gamma_high", 47, "-3.92", "bf10", "86.76"),
    ("P2 PLV ROI gamma_high", 47, "-3.92", "dz", "-0.57"),
    ("P2 ImC ROI", 47, "-0.146", "p", "0.884"),
    ("P2 ImC ROI", 47, "-0.146", "bf01", "6.25"),
    ("P2 ImC ROI", 47, "-0.146", "p_tost", "0.0063"),
    ("P2 ImC F2-P4", 47, "0.259", "p", "0.797"),
    ("P2 ImC F2-P4", 47, "0.259", "bf01", "6.12"),
    ("P2 ImC F2-P4", 47, "0.259", "p_tost", "0.0084"),
    ("P2 CSD ROI", 47, "-0.63", "bf01", "5.23"),
    ("P2 CSD ROI", 47, "-0.63", "p_tost", "0.020"),
    ("P2 CSD F2-P4", 47, "1.39", "bf01", "2.58"),
    ("P2 ImC under CSD", 47, "-2.54", "bf10", "2.79"),
    ("ceiling N=47", 47, "0", "bf01", "6.32"),   # handoff said 6.31: that is r = 0.707,
                                                  # inconsistent with its own 4.49 at N = 22
    ("ceiling N=22", 22, "0", "bf01", "4.49"),
    ("min TOST p N=22", 22, "0", "p_tost", "0.037"),
]


def _half_ulp(s: str) -> float:
    return 0.5 * 10.0 ** -(len(s.split(".")[1]) if "." in s else 0) if s != "0" else 0.0


def _selftest() -> bool:
    ok_all = True
    print(f"{'check':24s} {'field':7s} {'printed':>8s} {'computed range':>21s}  result")
    for label, n, ts, field, rep in _TARGETS:
        t, ht = float(ts), _half_ulp(ts)
        vals = [summary_from_t(t + d, n)[field] for d in (-ht, 0.0, ht)]
        lo, hi = min(vals), max(vals)
        r, hr = float(rep), _half_ulp(rep)
        ok = (hi >= r - hr - 1e-12) and (lo < r + hr)
        ok_all &= ok
        print(f"{label:24s} {field:7s} {rep:>8s} [{lo:9.5f}, {hi:9.5f}]  {'PASS' if ok else 'FAIL'}")

    rng = np.random.default_rng(0)
    x = rng.normal(0.3, 1.0, 10)
    sx2, n = float(np.sum(x * x)), x.size
    t_obs = abs(float(_t_from_sum(x.sum(), sx2, n)))
    brute = np.mean([abs(float(_t_from_sum(np.dot(s, x), sx2, n))) >= t_obs - 1e-10
                     for s in itertools.product((-1.0, 1.0), repeat=n)])
    p_exact, _ = signflip_p(x)
    ok = abs(brute - p_exact) < 1e-12
    ok_all &= ok
    print(f"{'sign-flip exact vs brute':24s} {'p_perm':7s} {brute:8.5f} {p_exact:21.5f}  {'PASS' if ok else 'FAIL'}")
    return ok_all


if __name__ == "__main__":
    raise SystemExit(0 if _selftest() else 1)
