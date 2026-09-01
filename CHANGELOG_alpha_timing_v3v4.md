# Alpha-validation timing: v3 → v4 correction

Referenced from `config.m` (`cfg.alpha.lock_event = 'cue'`). This document is
historical record of an already-applied decision, not a pending patch — all
the parameter values described below are already live in `config.m`
(`cfg.ti.col_cue_sample`, `cfg.alpha.*`).

## v4 correction (supersedes v3)

The original paper (Gregory, Wang & Kessler, 2022, *SCAN*) states explicitly
in Methods: *"epoched from 1s pre cue onset to 1s post probe response, such
that CUE ONSET = TIME 0."* Figure 1's timeline confirms: Targets shown
3.5–4s, Maintenance 4–5s, Probe 5–6s — **all relative to cue onset**, not
trial Start.

v3's empirical timing check correctly identified the column mapping
(Maint → LocProbe ≈ 1s, matching the paper's Maintenance interval) but
**wrongly concluded t=0 was trial Start**. It is actually Cue onset — Start
is simply the ~1s fixation/ITI period before the cue (matches the paper's
stated 1000ms ITI). Locking to Start introduced a constant ~1.04s offset
into every analysis window, which is why the v3 run found:

- **(a)** no effect in what was meant to be the encoding window — it was
  actually analyzing the 2–3s post-cue "eye contact hold" period, *before*
  cue validity is even revealed to the participant at 3–3.5s, and
- **(b)** a spurious "retrieval" effect that was actually sitting right on
  the target-offset/maintenance-onset transition (~cue+4.0s) — an evoked-response
  artifact, not a retrieval-interval effect.

(The diagnostic scripts that surfaced (b) — `s10_diagnose_retrieval_effect.m`
and `s11_diagnose_encoding_topo.m` — are kept in `/alpha_validation/diagnostics/`
as the paper trail for this finding, not as part of the live pipeline.)

## What v4 changed

v4 reverts to locking on cue onset directly (`cfg.ti.col_cue_sample`).
Because the epoch's own t=0 *is* each trial's actual cue-onset sample,
cue-onset jitter is architecturally irrelevant here (unlike the v3
Start-locked design) — so the per-trial manual baseline correction added in
v3.1 is no longer needed either; a single fixed baseline window (matching
the paper's own −500 to −100ms pre-cue) is now exactly correct, not an
approximation.

## Where this lives now

All parameters from this investigation (`cfg.ti.col_start/col_cue_sample/col_targets_sample`,
`cfg.alpha.lock_event/baseline_s/analysis_s/foi/cycles/toi_step_s/band/avgoverfreq/n_perm/stat_latency_*`)
are defined once, in `config.m`. This file exists only to explain *why*
those specific values were chosen.
