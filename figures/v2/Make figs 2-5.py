"""Regenerate Figures 2-5 of the Pilot 2 preprint (ds003702) from results/.

Usage:  python make_figs_2_5.py <path-to-repo>/results <output-dir>

Every plotted value is read from the results/ CSVs of the post-Muscle-IC-fix
run; nothing is typed in by hand. The script asserts the headline values
against the preprint text, so a stale or mismatched results folder fails
loudly instead of producing a plausible-looking figure.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import integrate

RES = Path(sys.argv[1])
OUT = Path(sys.argv[2])
OUT.mkdir(parents=True, exist_ok=True)

# Palette: validated (dataviz validate_palette.js, light mode, all checks pass)
C_EFFECT = "#eb6834"   # PLV: evidence for an effect
C_NULL = "#2a78d6"     # robust measures: evidence for the null
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8985"
GRID = "#e4e3df"
SURF = "#ffffff"

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2,
    "xtick.color": INK2, "ytick.color": INK2,
    "axes.spines.top": False, "axes.spines.right": False,
    "savefig.facecolor": SURF, "figure.facecolor": SURF,
})


def jzs_bf10(t, n, r=0.707):
    """One-sample JZS Bayes factor (Rouder et al., 2009), two-sided."""
    nu = n - 1
    f = lambda g: ((1 + n * g * r**2) ** -0.5
                   * (1 + t**2 / ((1 + n * g * r**2) * nu)) ** (-(nu + 1) / 2)
                   * (2 * np.pi) ** -0.5 * g ** -1.5 * np.exp(-1 / (2 * g)))
    num, _ = integrate.quad(f, 0, np.inf, limit=500)
    return num / (1 + t**2 / nu) ** (-(nu + 1) / 2)


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"{name}.{ext}", dpi=300, bbox_inches="tight")
    plt.close(fig)


# ----------------------------------------------------------------- load
bt = pd.read_csv(RES / "bayes_tost_gamma.csv")
row = lambda m: bt[(bt.contrast == "d_val") & (bt.metric == m)
                   & (bt.band == "gamma_high")].iloc[0]
plv, dwpli = row("plv"), row("dwpli")
imc_roi = pd.read_csv(RES / "imag_coherence_roi_average_gamma_high.csv").iloc[0]
imc_pairs = pd.read_csv(RES / "imag_coherence_valid_main_gamma_high.csv")
imc_f2p4 = imc_pairs[imc_pairs.pair == "F2-P4"].iloc[0]

sens = pd.read_csv(RES / "sensitivity_gammahigh_bayes.csv")
s_sub48 = sens[(sens.policy == "WITHOUT sub-48") & (sens.contrast == "d_val")
               & (sens.metric == "plv") & (sens.band == "gamma_high")].iloc[0]
rts = pd.read_csv(RES / "sensitivity_rtstat_gammahigh.csv")
s_rt = rts[rts.policy.str.startswith("RT_stat filter DISABLED")
           & (rts.contrast == "d_val") & (rts.metric == "plv")].iloc[0]

ph = pd.read_csv(RES / "phase_diff_diagnostic_F2P4_gamma_high.csv")

# Guard against a stale results folder: values must match the preprint text.
assert round(plv.BF10, 2) == 86.76, plv.BF10
assert round(dwpli.BF01, 2) == 6.18, dwpli.BF01
assert round(imc_roi.BF01, 2) == 6.25, imc_roi.BF01
assert round(imc_f2p4.BF01, 2) == 6.12, imc_f2p4.BF01
assert round(s_sub48.BF10, 2) == 86.91 and round(s_rt.BF10, 2) == 87.51
assert len(ph) == 188

N = int(plv.n)
bf01_ceiling = 1 / jzs_bf10(0.0, N)     # max attainable BF01 at this N (t = 0)

# ------------------------------------------------------------ Figure 3
rows = [
    ("PLV\nROI average", plv.BF10, f"BF$_{{10}}$ = {plv.BF10:.2f}", C_EFFECT),
    ("dwPLI\nROI average", 1 / dwpli.BF01, f"BF$_{{01}}$ = {dwpli.BF01:.2f}", C_NULL),
    ("Imaginary coherency\nROI average", 1 / imc_roi.BF01,
     f"BF$_{{01}}$ = {imc_roi.BF01:.2f}", C_NULL),
    ("Imaginary coherency\nF2–P4", 1 / imc_f2p4.BF01,
     f"BF$_{{01}}$ = {imc_f2p4.BF01:.2f}", C_NULL),
]
fig, ax = plt.subplots(figsize=(6.5, 2.9))
ax.set_xscale("log")
ax.set_xlim(0.03, 400)
ax.axvspan(1 / 3, 3, color="#f1f0ec", zorder=0, lw=0)
ax.axvline(1, color=MUTED, lw=1, zorder=1)
ax.axvline(1 / bf01_ceiling, color=INK2, lw=1, ls=(0, (2, 2)), zorder=1)
ys = np.arange(len(rows))[::-1]
# Bars start at BF = 1 (no evidence either way) and extend left or right.
for y, (lab, bf10, txt, col) in zip(ys, rows):
    ax.barh(y, bf10 - 1, left=1, height=0.55, color=col, zorder=2)
    if bf10 > 1:
        ax.text(bf10 * 1.12, y, txt, va="center", ha="left", color=INK, fontsize=8.5)
    else:  # label inside the bar, white on blue (contrast > 4.5:1)
        ax.text(0.92, y, txt, va="center", ha="right", color="white", fontsize=8.5)
ax.set_yticks(ys, [r[0] for r in rows], fontsize=8.5)
ax.tick_params(axis="y", length=0)
ticks = [0.1, 1 / 3, 1, 3, 10, 30, 100]
ax.set_xticks(ticks, ["1/10", "1/3", "1", "3", "10", "30", "100"])
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
ax.set_xlabel("Bayes factor BF$_{10}$ (log scale)   ← evidence for null  |  evidence for effect →")
ax.text(1 / bf01_ceiling / 1.08, ys[0] + 0.45,
        f"max. attainable\nBF$_{{01}}$ at N = {N}\n({bf01_ceiling:.2f})",
        fontsize=7.5, color=INK2, va="bottom", ha="right")
ax.text(1, ys[0] + 0.45, "inconclusive\n(1/3 to 3)", fontsize=7.5, color=INK2,
        va="bottom", ha="center")
ax.set_ylim(-0.6, len(rows) + 0.35)
save(fig, "Figure3_bayes_factors")

# ------------------------------------------------------------ Figure 4
ang = np.deg2rad(ph.mean_angle_deg.values)
rad = ph.resultant_len.values
dev = np.abs((ph.mean_angle_deg.values % 360) - 180)
inside = dev <= 20
n_in, n_tot = int(inside.sum()), len(ph)
assert n_in == 175, n_in

fig = plt.figure(figsize=(4.6, 4.6))
ax = fig.add_subplot(projection="polar")
ax.set_facecolor(SURF)
ax.bar(np.pi, 0.75, width=np.deg2rad(40),
       color="#f6e7b0", alpha=0.6, zorder=0, lw=0)
ax.scatter(ang[inside], rad[inside], s=16, color=C_NULL, edgecolor=SURF,
           lw=0.6, zorder=3, label=f"within ±20° of 180° (n = {n_in})")
ax.scatter(ang[~inside], rad[~inside], s=16, color=MUTED, edgecolor=SURF,
           lw=0.6, zorder=4, label=f"outside ±20° of 180° (n = {n_tot - n_in})")
ax.set_ylim(0, 0.75)
ax.set_yticks([0.2, 0.4, 0.6])
ax.set_yticklabels(["0.2", "0.4", "0.6"], fontsize=7.5, color=INK2)
ax.set_rlabel_position(60)
ax.set_xticks(np.deg2rad([0, 45, 90, 135, 180, 225, 270, 315]))
ax.grid(color=GRID, lw=0.6)
ax.spines["polar"].set_color(MUTED)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.06), ncol=1,
          frameon=False, fontsize=8, handletextpad=0.3)
save(fig, "Figure4_phase_difference_F2P4")

# ------------------------------------------------------------ Figure 5
srows = [
    ("Primary analysis\n(N = 47)", plv),
    ("Excluding sub-48\n(N = 46)", s_sub48),
    ("RT_stat filter disabled\nfor sub-01 to sub-04", s_rt),
]
fig, ax = plt.subplots(figsize=(6.5, 2.2))
ax.set_xscale("log")
ax.set_xlim(1, 300)
ax.axvspan(1, 3, color="#f1f0ec", zorder=0, lw=0)
ys = np.arange(len(srows))[::-1]
for y, (lab, r) in zip(ys, srows):
    ax.hlines(y, 1, r.BF10, color=C_EFFECT, lw=2, zorder=2)
    ax.scatter(r.BF10, y, s=46, color=C_EFFECT, edgecolor=SURF, lw=1.2, zorder=3)
    ax.text(r.BF10 * 1.15, y,
            (f"BF$_{{10}}$ = {r.BF10:.2f}   (d$_z$ = {r.dz:.2f}, "
             f"t({int(r.n) - 1}) = {r.t:.2f})").replace("-", "\u2212"),
            va="center", ha="left", fontsize=8.5, color=INK)
ax.set_yticks(ys, [s[0] for s in srows], fontsize=8.5)
ax.tick_params(axis="y", length=0)
ax.set_xticks([1, 3, 10, 30, 100], ["1", "3", "10", "30", "100"])
ax.xaxis.set_minor_locator(matplotlib.ticker.NullLocator())
ax.grid(axis="x", color=GRID, lw=0.6, zorder=0)
ax.set_xlabel("Bayes factor BF$_{10}$, validity effect on high-gamma PLV (log scale)")
ax.set_ylim(-0.6, len(srows) - 0.4)
save(fig, "Figure5_sensitivity")

print(f"BF01 ceiling at N={N}: {bf01_ceiling:.3f}")
print(f"Phase cells within 20 deg of 180: {n_in}/{n_tot} = {100 * n_in / n_tot:.1f}%")
print("PLV", plv.BF10, "dwPLI", dwpli.BF01, "ImC ROI", imc_roi.BF01, "ImC F2-P4", imc_f2p4.BF01)
print("TOST p: dwPLI", dwpli.p_TOST, "ImC ROI", imc_roi.p_TOST, "ImC F2-P4", imc_f2p4.p_TOST)
print("Sensitivity:", s_sub48.BF10, s_rt.BF10)


# ------------------------------------------------------------ Figure 2
# Frontoparietal validity cluster (high-gamma PLV), from s06/s06c outputs.
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Circle, Polygon, Rectangle, Ellipse
from matplotlib.path import Path as MPath
from matplotlib.patches import PathPatch

members = pd.read_csv(RES / "cluster_mass_valid_main_plv_gamma_high_members.csv")
cmass = pd.read_csv(RES / "cluster_mass_valid_main_plv_gamma_high.csv")
pairs = pd.read_csv(RES / "conn_cluster_results_valid_main_by_metric.csv")
pairs = pairs[(pairs.metric == "plv") & pairs.pair_band.str.startswith("gamma_high")].copy()
pairs["pair"] = pairs.pair_band.str.replace("gamma_high ", "", regex=False)
tmap = dict(zip(pairs.pair, pairs.t)); pfw = dict(zip(pairs.pair, pairs.p_fwer))

clus = members[members.cluster_id == 1]
c1 = cmass[cmass.cluster_id == 1].iloc[0]
strongest = pairs.loc[pairs.t.abs().idxmax(), "pair"]
diag_pair = "F2-P4"
assert strongest == "FCz-Pz" and len(clus) == 29 and round(c1.p_cluster, 4) == 0.0002
assert diag_pair in set(clus.pair)

FRONT = ["F3", "F1", "Fz", "F2", "F4", "FC1", "FCz", "FC2"]   # cfg.roi_frontal (display order)
PAR = ["CP1", "CP2", "P3", "P1", "Pz", "P2", "P4", "POz"]      # cfg.roi_parietal (display order)
assert len(pairs) == 64
# Approximate 2-D scalp positions (polar angle from nose in deg, radius with
# head edge = 0.5), standard 10-10 layout; schematic only.
POS = {"F3": (-40, .34), "F1": (-22, .27), "Fz": (0, .25), "F2": (22, .27), "F4": (40, .34),
       "FC1": (-42, .17), "FCz": (0, .13), "FC2": (42, .17),
       "CP1": (-138, .17), "CP2": (138, .17), "P3": (-140, .34), "P1": (-158, .27),
       "Pz": (180, .25), "P2": (158, .27), "P4": (140, .34), "POz": (180, .38)}
XY = {k: (r * np.sin(np.deg2rad(a)), r * np.cos(np.deg2rad(a))) for k, (a, r) in POS.items()}
# The four midline-to-midline pairs are collinear; give each a distinct bow.
BOW = {"Fz-Pz": .05, "Fz-POz": -.05, "FCz-POz": .085, "FCz-Pz": 0.0}

div = LinearSegmentedColormap.from_list(
    "div", ["#104281", "#2a78d6", "#f0efec", "#e34948", "#8f1f1f"])
TMAX = 4.0

fig = plt.figure(figsize=(7.0, 3.6))
axh = fig.add_axes([0.0, 0.12, 0.42, 0.84])
axm = fig.add_axes([0.53, 0.20, 0.36, 0.70])
cax = fig.add_axes([0.905, 0.20, 0.015, 0.70])

# (a) head schematic
axh.set_aspect("equal"); axh.axis("off")
axh.add_patch(Circle((0, 0), 0.5, fill=False, ec=MUTED, lw=1.2))
axh.add_patch(Polygon([(-0.05, 0.495), (0, 0.56), (0.05, 0.495)], closed=False,
                      fill=False, ec=MUTED, lw=1.2))
for sx in (-1, 1):
    axh.add_patch(Ellipse((sx * 0.515, 0), 0.04, 0.14, fill=False, ec=MUTED, lw=1.2))


def arc(a, b, bow, **kw):
    (x1, y1), (x2, y2) = XY[a], XY[b]
    mx, my = (x1 + x2) / 2, (y1 + y2) / 2
    dx, dy = x2 - x1, y2 - y1; L = np.hypot(dx, dy)
    cx, cy = mx - bow * dy / L, my + bow * dx / L
    path = MPath([(x1, y1), (cx, cy), (x2, y2)], [MPath.MOVETO, MPath.CURVE3, MPath.CURVE3])
    axh.add_patch(PathPatch(path, fill=False, capstyle="round", **kw))


for pr in clus.pair:
    if pr in (strongest, diag_pair):
        continue
    f, p_ = pr.split("-")
    arc(f, p_, BOW.get(pr, 0.0), ec="#b9b8b2", lw=0.9, zorder=1)
f, p_ = diag_pair.split("-"); arc(f, p_, 0.0, ec=INK, lw=2.0, ls=(0, (3, 2)), zorder=3)
f, p_ = strongest.split("-"); arc(f, p_, BOW.get(strongest, 0.0), ec=INK, lw=2.4, zorder=3)

in_clus = {e for pr in clus.pair for e in pr.split("-")}
for e, (x, y) in XY.items():
    filled = e in in_clus
    axh.add_patch(Circle((x, y), 0.021, fc=INK2 if filled else SURF, ec=INK2, lw=1.0, zorder=4))
    dy = 0.042 if POS[e][1] > 0 and abs(POS[e][0]) < 90 else -0.045
    axh.text(x, y + dy, e, ha="center", va="center", fontsize=7, color=INK, zorder=5,
             bbox=dict(boxstyle="round,pad=0.12", fc=SURF, ec="none", alpha=0.9))
axh.set_xlim(-0.58, 0.58); axh.set_ylim(-0.56, 0.6)
axh.text(-0.58, 0.6, "a", fontsize=11, fontweight="bold", va="top")

# (b) t-matrix, all 64 frontal x parietal pairs
M = np.array([[tmap[f"{f}-{p}"] for p in PAR] for f in FRONT])
im = axm.imshow(M, cmap=div, vmin=-TMAX, vmax=TMAX, aspect="equal")
clus_set = set(clus.pair)
for i, f in enumerate(FRONT):
    for j, p in enumerate(PAR):
        pr = f"{f}-{p}"
        # 2px surface gap between cells
        axm.add_patch(Rectangle((j - .5, i - .5), 1, 1, fill=False, ec=SURF, lw=1.5, zorder=2))
        if pr in clus_set:
            axm.add_patch(Rectangle((j - .42, i - .42), .84, .84, fill=False, ec=INK2, lw=0.8, zorder=3))
        if pr == strongest:
            axm.add_patch(Rectangle((j - .45, i - .45), .9, .9, fill=False, ec=INK, lw=2.2, zorder=4))
        if pr == diag_pair:
            axm.add_patch(Rectangle((j - .45, i - .45), .9, .9, fill=False, ec=INK, lw=1.8,
                                    ls=(0, (3, 2)), zorder=4))
axm.set_xticks(range(len(PAR)), PAR, fontsize=7.5, rotation=0)
axm.set_yticks(range(len(FRONT)), FRONT, fontsize=7.5)
axm.tick_params(length=0)
for s_ in axm.spines.values():
    s_.set_visible(False)
axm.set_xlabel("Parietal", fontsize=8.5); axm.set_ylabel("Frontal", fontsize=8.5)
axm.text(-0.2, 1.08, "b", transform=axm.transAxes, fontsize=11, fontweight="bold", va="top")
cb = fig.colorbar(im, cax=cax)
cb.outline.set_visible(False); cb.ax.tick_params(labelsize=7, length=2)
cb.set_label("t(46), valid − invalid", fontsize=8)
cb.ax.yaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}".replace("-", "−")))

handles = [
    Line2D([], [], color=INK, lw=2.4,
           label=f"{strongest.replace('-', '–')}: strongest pair (p$_{{FWER}}$ = {pfw[strongest]:.3f})"),
    Line2D([], [], color=INK, lw=1.8, ls=(0, (3, 2)),
           label=f"{diag_pair.replace('-', '–')}: diagnostic pair (p$_{{FWER}}$ = {pfw[diag_pair]:.3f})"),
    Line2D([], [], color="#b9b8b2", lw=1.2,
           label=f"other members of the {len(clus)}-pair cluster (p$_{{cluster}}$ = {c1.p_cluster:.4f})"),
]
fig.legend(handles=handles, loc="lower center", ncol=1, frameon=False, fontsize=7.5,
           bbox_to_anchor=(0.5, -0.04))
save(fig, "Figure2_validity_cluster")
print("Fig2:", strongest, pfw[strongest], diag_pair, pfw[diag_pair], len(clus), c1.p_cluster,
      "electrodes not in cluster:", sorted(set(XY) - in_clus))
