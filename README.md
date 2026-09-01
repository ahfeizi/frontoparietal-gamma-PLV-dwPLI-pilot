# A PLV–dwPLI Dissociation in Gamma-Band Frontoparietal Connectivity — ds003702 (Social Memory Cuing)

A methodological reanalysis (N = 47) of a public VR social-cueing working-memory
dataset (OpenNeuro ds003702; Gregory, Wang & Kessler, 2022), testing whether
frontoparietal **gamma-band** phase synchronization is modulated by cue-target
validity — and whether an apparent PLV effect survives volume-conduction-robust
scrutiny (dwPLI, imaginary coherency, trial-level phase-difference distributions).

**Headline finding:** a strong, spatially coherent PLV effect (BF10 = 86.76,
p_cluster = 0.0002 across 29 pairs) does **not** survive any
volume-conduction-robust check (dwPLI BF01 = 6.18; imaginary coherency
BF01 = 6.12–6.25; 93% of phase differences within 20° of 180°; collapses
under CSD re-referencing to BF01 = 5.23 ROI-average). This dissociation — not
the PLV effect itself — is the paper's contribution, and directly informs the
volume-conduction-aware connectivity strategy adopted for the planned
tACS-EEG-fMRI study.

**Post-fix note:** all numbers on this page reflect the corrected
`s01_preprocess_New.m` pipeline (dominant-channel Muscle-IC interpolation
bug fixed — see `/archive/` note below) and are cross-checked against the
revised preprint. The PLV effect got *stronger* after the fix (BF10:
45.16→86.76), not weaker — see the verification log in this repo's PR
history / commit messages for the full pre-fix vs. post-fix comparison.

This pilot is the companion analysis to
[`ahfeizi/frontotemporal-WM-pilot`](https://github.com/ahfeizi/frontotemporal-WM-pilot)
(Pilot 1, ds004117) and follows the same open-science / verification conventions.

---

## Repository structure

```
/
  config.m                          <- central config; EDIT paths before running
  CHANGELOG_alpha_timing_v3v4.md    <- why alpha-validation locks to cue onset, not trial start
  README.md
  /pipeline/                        <- core connectivity pipeline (produces the headline result)
  /behavioral/                      <- behavioral replication (Table 1)
  /alpha_validation/                <- independent replication of the original paper's alpha effect
    /diagnostics/                   <- QC/exploratory scripts behind that replication, not load-bearing
  /connectivity_gamma/               <- everything downstream of s03: stats, Bayes/TOST, dissociation,
                                        sensitivity, CSD checks
  /figures/
  /results/                         <- numeric outputs cited directly in the preprint (see table below)
  /archive/                         <- superseded script versions, kept for provenance
```

## Run order

```
/pipeline/
  s00_inventory              % fast, single-subject channel/event check (v2: no full pop_importbids)
  s01_preprocess_New         % per-subject preprocessing + manual bad-channel override (sub-48/P1)
  s02_epoch                  % trial list from the AUTHORS' trialinfo (not re-parsed from raw markers)
  s03_connectivity           % THE SCIENTIFIC CORE: dwPLI + PLV, gamma_low/gamma_high, all 2x2 cells

/behavioral/
  s00b_phase0_behavior       % PRIMARY: repeated-measures ANOVA, 4 DVs x 3 effects -> Table 1 source
  s07_behavioral              % independent cross-check via permutation test (subset of the same effects)

/alpha_validation/
  s08_alpha_epoch_tfr        % per-subject cue-locked alpha TFR (v4: locked to cue onset, not trial start)
  s08_patch_add_elec         % one-time patch (only needed if s08_alpha_epoch_tfr ran pre-elec-fix)
  s09_alpha_stats            % group-level cluster permutation, restricted to the paper's own windows
  /diagnostics/
    s02b_check_maint_timing, s08a_check_trial_timeline, s08a_explore_cue_markers,
    s10_diagnose_retrieval_effect (rename: diag_retrieval_effect),
    s11_diagnose_encoding_topo (rename: diag_encoding_topo),
    s12_evoked_induced_decomposition (rename: diag_evoked_induced)

/connectivity_gamma/
  s04_stats                   % cue_main / valid_main / interaction, sign-flip permutation
  s06_cluster_conn             % FWER max-stat test, full pair space, dwPLI + PLV, 3 contrasts
  s06c_robustness_valid_plv_gammahigh   % 5 seeds x {5000,10000} perms stability check on the F2-P4 cluster
  s08_trial_counts             % per-subject/cell trial counts (SNR/dwPLI-bias relevant)
  s10_bayes_tost_gamma         % BF10/BF01 (JZS) + TOST, all 3 contrasts x 2 metrics x 2 sub-bands
  s13_imaginary_coherence      % dissociation test: ImC vs PLV/dwPLI (+ phase_diff_diagnostic_F2P4)
  s14_sensitivity_sub48        % with/without sub-48 -- see "Known reproducibility gaps" below
  s15_sensitivity_rtstat       % RT_stat filter reactivated for sub-01..04
  s16_multiple_comparisons_audit  % Holm-Bonferroni, Family 1 (behavioral) + Family 2 (cluster)
  s17_visual_confound_phase4   % occipital ERP, avatar vs stick, 70-340ms
  s18_csd_sensitivity          % CSD re-reference, primary PLV/dwPLI/ImC contrast
  s19_csd_perpair_imc          % per-pair follow-up on the CSD ImC anomaly
  s20_csd_imc_cluster_followup % seed-stability (10 combos) + cross-metric check on that 7-pair cluster

/figures/
  s05_figures                  % ⚠️ predates the dissociation finding -- see note below
```

## `s05_figures.m` — rewritten

The original version (now in `/archive/`) plotted PLV and dwPLI as separate,
independent paired comparisons — predating the dissociation finding. The
rewritten version reads directly from verified `/results/` files (not
`conn_results.mat`) and produces:

- **Fig. 3** — BF10/BF01 across PLV, dwPLI, and imaginary coherency (the
  paper's central "multi-metric convergence" figure)
- **Fig. 4** — circular distribution of F2-P4 phase differences (188 cells)
- **Fig. 5** — stability of the PLV BF10 across the primary analysis and both
  sensitivity checks

**One open item inside the new script:** Fig. 3's imaginary-coherency bar
currently uses the F2-P4-specific value (BF01=6.31) because the ROI-averaged
ImC value the preprint text cites (BF01=6.18) isn't in any `results/` file
yet — only the per-pair `imag_coherence_valid_main_gamma_high.csv` exists.
The script flags this at runtime (`fprintf` warning) and in an inline
comment. If an ROI-averaged ImC summary file exists (or can be produced,
analogous to how `csd_sensitivity_summary.csv` has both `roi_avg` and
`F2P4` rows), send it and I'll wire it in.

**Not generated by this script:** Fig. 1 (behavioral replication — needs the
*original study's* F-statistics as reference points, which are only in the
preprint text right now, not in any results file) and Fig. 2 (cluster
schematic — a hand-drawn/Visualizer network diagram, not a data plot).

## Reproducibility gaps — status

All previously-flagged gaps are now resolved:

- **`s14`/`s15` (BF10: 45.16→45.09 and 45.16→46.05):** saved to
  `sensitivity_gammahigh_bayes.csv`, `sensitivity_f2p4_cluster.csv`, and
  `sensitivity_rtstat_gammahigh.csv`; verified against the preprint.
- **Visual confound cluster stat:** `visual_confound_timecourse.mat`
  (`Ma`, `Ms`: per-subject avatar/stick time courses; `tobs_tc`: observed
  t-time-course; `srate=500`) is now in `results/`. A direct check of the
  70–340ms window gives max|t|=10.03 within-window and a window-mean paired
  t-test of t(46)=7.82, p≈5×10⁻¹⁰ — this is a simple paired t on the window
  mean, not the actual cluster-permutation statistic the preprint reports,
  but it strongly corroborates the direction and magnitude of "p < 0.001."
- **`s05_figures.m`:** rewritten (see below).

## `results/` — numeric outputs cited directly in the preprint

All files below have been cross-checked against the **post-fix, final**
preprint numbers and match. (An earlier version of this pipeline had a bug
in `s01_preprocess_New.m` where dominant-channel Muscle-IC components were
never actually interpolated out. 37 of 47 subjects were affected. All
numbers here reflect the corrected pipeline, rerun end to end.)

| File | Produced by | Final verified claim |
|---|---|---|
| `phase_diff_diagnostic_F2P4_gamma_high.csv` | `s13` | Fig. 4: 175/188 cells within 20° of 180° (93%) |
| `bayes_tost_gamma.csv` | `s10` | BF10=86.76 (PLV), BF01=6.18 (dwPLI) |
| `imag_coherence_valid_main_gamma_high.csv` | `s13` | F2-P4: BF01=6.12, dz=0.038 |
| `imag_coherence_roi_average_gamma_high.csv` | `s13` (patched — see below) | ROI-avg: BF01=6.25 |
| `cluster_mass_valid_main_plv_gamma_high.csv` + `_members.csv` | `s06c` | 29-pair cluster, p_cluster=0.0002 (FieldTrip montage adjacency — the only adjacency definition retained; see limitations below) |
| `conn_cluster_results_{cue_main,valid_main,interaction}_by_metric.csv` | `s06` | full pairwise max-stat space; strongest single pair is now FCz-Pz (p_fwer=0.0512), not F2-P4 (p_fwer=0.0846) |
| `multcomp_audit_family1_behavior.csv` | `s16` | Table 1 — unchanged by the s01 fix (behavioral data, as expected) |
| `multcomp_audit_family2_cluster.csv` | `s16` | strongest pair FCz-Pz, p_fwer=0.0512, p_holm=0.3072 |
| `csd_sensitivity_summary.csv` + `_persubject.csv` | `s18` | ROI-avg: BF10 86.76→BF01=5.23; F2-P4 t: -3.52→+1.39 (BF10=2.58, not equivalent); ROI-avg ImC under CSD: BF01=2.79 (inconclusive, not a positive signal) |
| `csd_imc_perpair_valid_main_gamma_high.csv` | `s19` | per-pair breakdown underlying the CSD ImC cluster |
| `csd_imc_clustermass_valid_main_gamma_high.csv` | `s19` | **7-pair** cluster (not 9), mass=17.12, p_cluster=0.023 |
| `csd_imc_cluster_seedstability.csv` | `s20` (patched — see below) | stable across 10 seed/nperm combos, p_cluster 0.023–0.027 |
| `csd_imc_cluster_crossmetric.csv` | `s20` (patched) | 0 of 7 pairs significant in dwPLI/PLV (not corroborated) |
| `sensitivity_acc_loc_FULL__all_47_.csv` / `_WITHOUT_sub-48.csv` | `s14` (Part C) | behavioral robustness to sub-48 exclusion — unchanged by s01 fix |
| `sensitivity_gammahigh_bayes.csv` | `s14` (Parts A/B, patched) | BF10 with/without sub-48 |
| `sensitivity_f2p4_cluster.csv` | `s14` (patched) | matches `s06c`'s 29-pair cluster exactly (78.51, p=0.0002) — cross-script consistency check passed |
| `sensitivity_rtstat_gammahigh.csv` | `s15` (patched) | BF10 with/without RT_stat reactivation |
| `visual_confound_occipital_rms.csv` + `visual_confound_timecourse.mat` | `s17` | descriptive RMS + time-resolved t-course; window-mean sanity check t(46)=7.82, p≈5e-10, consistent with the reported p<0.001 cluster stat |
| `conn_results_rtstat_sensitivity.mat` | `s15` | raw recomputed connectivity for sub-01..04 |

**Two cache/staleness bugs were found and fixed during verification** (see
commit history for exact diffs):
1. `s01_preprocess_New.m` — dominant-channel Muscle-IC interpolation never
   actually executed (`interp_channels` was computed but never fed into a
   `pop_select` before the final `pop_interp`). Fixed; affected 37/47
   subjects; the full pipeline (`s02`→`s20`) was rerun from the corrected
   preprocessing output.
2. `s20_csd_imc_cluster_followup.m` — silently reused a stale cached D
   matrix (`csd_imc_dmatrix_cache.mat`) and hardcoded pre-fix cluster
   membership/mass, producing numbers from the *old* preprocessing even
   after `s19` had already been correctly rerun. Fixed (stale cache
   deleted, hardcoded values updated); the corrected 7-pair/mass=17.12
   result now matches `s19`'s independent recomputation exactly.

Three additional bugs were found and fixed but did not change any
downstream numbers: `config.m` hardcoded personal paths; a personal
`addpath` in `s03_connectivity.m` exposing the host supervisor's name; and
`s06c_robustness_valid_plv_gammahigh.m`'s `neighbours_file` being assigned
before `cfg = config()` existed (would error on a fresh clone).

## `/archive/`

- `s01_preprocess.m` — superseded by `s01_preprocess_New.m`, which adds the
  `cfg.manual_bad_channels` override (sub-48/P1 interpolation) that the rest
  of the pipeline assumes is applied. Running the archived version would
  silently invalidate the sub-48 sensitivity analysis (`s14`).
- `config_old.m` — superseded by `config.m`.
- `s08b_alpha_validation_pipeline.m` — an earlier, monolithic attempt at the
  alpha-validation replication; never run to completion. Superseded by the
  split `s08_alpha_epoch_tfr.m` + `s09_alpha_stats.m` pipeline.

## Honest limitations (from the preprint, stated plainly)

- This is a reanalysis of a single public dataset; generalizability of the
  specific PLV artifact to other montages, references, or gamma windows is
  not established.
- dwPLI and imaginary coherency cannot, by themselves, distinguish volume
  conduction from an average-reference artifact. CSD re-referencing (§3.9)
  rules out a purely reference-independent account but cannot on its own
  isolate which of the two zero-lag mechanisms is responsible.
- The exploratory 7-pair CSD/ImC cluster (`s19`/`s20`) is hypothesis-generating
  only: single-metric, present under only one reference scheme, identified
  post hoc, and not corrected for the additional comparisons this sensitivity
  branch entailed. Do not treat it as a confirmed finding.
- The cue × validity interaction in gamma connectivity remains genuinely
  inconclusive (BF01 in the 3.38–3.68 range, post-fix) — not a confirmed null.
- After Holm-Bonferroni correction across the six-test primary EEG family, no
  test survives, including the PLV effect (p_holm = 0.30). The paper's
  interpretive claim rests on multi-metric convergence, not a single
  corrected threshold — a deliberate departure from single-p-value reporting
  that should stay explicit wherever this pilot is cited (e.g. the MSCA
  proposal).

## Data and code availability

The preprint's current Data and Code Availability section has a placeholder:
*"[Repository link and DOI to be added upon deposition]"*. Once this repo is
public, update that sentence with the actual GitHub URL, and add the Zenodo
DOI once archived (matching the workflow already used for Pilot 1).

## Config

`config.m` requires local path edits (`cfg.bids_dir`, `cfg.eeglab_dir`,
`cfg.fieldtrip_dir`, and the `neighbours_file` path inside `s14`/`s06c`) before
running — none of these should point to a personal machine path once committed.
See `CHANGELOG_alpha_timing_v3v4.md` for the full investigation behind the
`cfg.alpha.lock_event = 'cue'` decision (§ referenced inline in `config.m`).
