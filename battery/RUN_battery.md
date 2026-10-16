# Running the convergent-validation battery (ds003702)

## Before you run

1. `config.m` must have a valid `cfg.deriv_dir`, `cfg.eeglab_dir` and `cfg.fieldtrip_dir`. These are existing fields; the **BATTERY ADDITIONS** block at the end adds only bands, wavelet frequencies, cycles and `cfg.battery.check_dir` (the folder the check scripts read from).
2. These existing outputs must already be in `cfg.deriv_dir`. Nothing is recomputed.
   - `manifest.mat`, `<key>_clean.set`/`.fdt` and `<key>_trials.mat` (from `s01_preprocess_New.m`, `s02_epoch.m`)
   - For the checks (`s22`, `s24`), read from `cfg.battery.check_dir`, which defaults to the repository's `results/` folder, where all of them are shipped. Set it to `cfg.deriv_dir` to check a fresh run instead:
     - `conn_results.mat` (`s03`)
     - `imag_coherence_roi_average_gamma_high.csv` (`s13`)
     - `csd_sensitivity_summary.csv` (`s18`)
     - `battery_persubject.csv`, `battery_phasediff.csv` (`s21`); `imc_diag_persubject.csv` (`s23`)
3. The scripts are in `battery/`. Each one adds the repository root and `battery/` to the MATLAB path, so `config.m`, `util_bandpass.m`, `util_chanidx.m` and `util_conn_pairs.m` all come from this repo. Set MATLAB's current folder to the repository root, because the relative paths in `config.m` (`cfg.deriv_dir`, `cfg.eeglab_dir`, `cfg.fieldtrip_dir`) are resolved from it.

## Toolboxes on the path

| Toolbox | Used for | How it is added |
|---|---|---|
| EEGLAB | `pop_loadset`, `eeglab2fieldtrip` | `addpath(cfg.eeglab_dir); eeglab nogui` (inside the scripts) |
| FieldTrip | `ft_scalpcurrentdensity` (CSD, method `'spline'`, same as `s18`) | `addpath(cfg.fieldtrip_dir); ft_defaults` (inside the scripts) |
| Signal Processing Toolbox | `fir1`, `filtfilt`, `hilbert` | MATLAB installation |

## Run order

1. `s21_battery.m` writes these files to `cfg.deriv_dir`:
   - `battery_persubject.csv`
   - `battery_phasediff.csv`
   - `battery_phasediff_bins.csv`
   - `battery_done_subjects.txt`
2. `s22_battery_check.m` prints PASS/FAIL for each check, then an `OVERALL` line.
3. `s23_imc_diagnostic.m` (exploratory, gamma_high only) writes `imc_diag_persubject.csv`, `imc_diag_power.csv` and `imc_diag_done_subjects.txt` to `cfg.deriv_dir`.
4. `s24_imc_diagnostic_check.m` checks `s23` against `battery_persubject.csv` and prints PASS/FAIL, then an `OVERALL` line.

## Resuming and restarting

`s21_battery.m` appends one subject at a time. After a crash, run it again:
- Subjects listed in `battery_done_subjects.txt` are skipped.
- Rows from a subject that was interrupted partway are removed before that subject is recomputed.

To start from scratch, delete the three `battery_*.csv` files and `battery_done_subjects.txt`.

## Expected runtime (estimate, not measured)

- **Per subject:** about 2–4 min (one continuous 500 Hz recording, 61 channels). The script prints the actual minutes for each subject.
- **Where the time goes:** loading the `.set`, the CSD on the full recording, and, for each reference, 21 Morlet bins on the 16 ROI channels.
- **All 47 subjects:** about 2–3 h.
- **Peak memory:** a few GB per subject, as in `s18` (continuous data plus its FieldTrip and CSD copies).
- **Check script:** a few minutes.
- **Output size:** `battery_persubject.csv` has about 578,000 rows (47 × 4 × 2 × 3 × 2 × 2 × 64).
