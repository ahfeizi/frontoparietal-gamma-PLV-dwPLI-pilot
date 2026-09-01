% s08_alpha_epoch_tfr.m  --  Per-subject cue-locked alpha-validation TFR
% (Pilot 2 / ds003702). Independent replication of the paper's cue-locked
% alpha effect, run BEFORE trusting the gamma dwPLI results from the
% maintenance-locked epochs (s02-s04).
%
% v4 CORRECTION (supersedes v3): the source paper (Gregory, Wang &
% Kessler, 2022, SCAN) states explicitly in Methods: "epoched from 1s
% pre cue onset to 1s post probe response, such that CUE ONSET = TIME 0."
% v3 had mistakenly locked to trial Start instead (based on an empirical
% timing check that identified the right columns but the wrong t=0
% reference). That silently shifted every analysis window by ~1s and
% produced (a) a null result in what was meant to be the encoding window
% -- it was actually analyzing the pre-cue-shift "eye contact hold"
% period, before validity is even revealed -- and (b) a spurious
% "retrieval" effect that was actually an evoked-response artifact at the
% target-offset/maintenance-onset transition. See config_ADDITIONS.m for
% the full explanation.
%
% This version locks to cue onset (col_cue_sample) directly, matching
% the paper exactly, and uses a single fixed baseline window (also
% matching the paper's own -500 to -100ms pre-cue exactly) -- the v3.1
% per-trial baseline machinery is no longer needed, since the epoch's
% own t=0 IS each trial's actual cue-onset sample, so jitter cannot
% misalign the baseline window regardless.
%
% ARCHITECTURE (matches config.m / s02_epoch.m): trial lock samples come
% directly from the AUTHORS' trialinfo absolute sample columns, NOT from
% re-parsing raw .vmrk marker pairs.
%
% For each subject:
%   1. Load our own cleaned continuous .set (s01_preprocess.m output) --
%      same file s02_epoch.m uses.
%   2. Load trialinfo from the authors' data_ica.mat (trialinfo only).
%   3. Recompute the SAME accuracy + RT inclusion filters as s02_epoch.m
%      (kept inline, independent of s02's saved trials.mat, since the
%      bounds/boundary check here uses a different window than the
%      maintenance-locked one).
%   4. Apply cfg.manual_bad_channels (same per-subject overrides).
%   5. Convert to FieldTrip, epoch on cfg.alpha.analysis_s relative to
%      cue onset, run Morlet TFR (2-30 Hz, 3 cycles, 50 ms steps),
%      dB-baseline correct (fixed window, exact per paper), split by
%      validity x cue-type, average.
%   6. Save <key>_alphaTFR.mat. Resume-capable via alpha_manifest.mat.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;

M = load(fullfile(cfg.deriv_dir, 'manifest.mat')); manifest = M.manifest;
srate = cfg.resample_hz;

win = round(cfg.alpha.analysis_s * srate);   % [start end] samples rel. to lock sample
toi = cfg.alpha.analysis_s(1):cfg.alpha.toi_step_s:cfg.alpha.analysis_s(2);
bl  = cfg.alpha.baseline_s;

if ~strcmp(cfg.alpha.lock_event, 'cue')
    error(['cfg.alpha.lock_event = ''%s'' but this script version assumes ''cue''. ' ...
        'If you intentionally changed it, update lock_col below accordingly.'], cfg.alpha.lock_event);
end
lock_col = cfg.ti.col_cue_sample;

alpha_manifest_file = fullfile(cfg.deriv_dir, 'alpha_manifest.mat');
done_keys = {};
if exist(alpha_manifest_file, 'file')
    dm = load(alpha_manifest_file); done_keys = dm.done_keys;
end

for r = 1:height(manifest)
    key  = char(manifest.key(r));
    subj = char(manifest.subject(r));
    out_file = fullfile(cfg.deriv_dir, [key '_alphaTFR.mat']);

    if ismember(key, done_keys) && exist(out_file, 'file')
        fprintf('%s already processed (alpha TFR), skipping (resume).\n', key);
        continue;
    end

    fprintf('\n=== %s (alpha validation, lock=%s) ===\n', key, cfg.alpha.lock_event);

    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);

    % NOTE: manual_bad_channels (e.g. P1 for sub-48) is NOT re-removed here.
    % s01_preprocess_New.m already excluded it before ICA and interpolated it
    % back at the end, so the loaded _clean.set already has a full, consistent
    % channel montage across all subjects -- removing it again here would
    % desync sub-48's channel count from every other subject going into
    % ft_freqstatistics (s09_alpha_stats.m uses subject 1's electrode layout
    % for everyone).

    ti_path = cfg.authors_deriv_file(subj);
    if ~exist(ti_path, 'file')
        warning('%s: trialinfo file not found at %s, skipping.', subj, ti_path);
        continue;
    end
    TI = load(ti_path, 'trialinfo');
    trialinfo = TI.trialinfo;
    nTrials = size(trialinfo, 1);

    lock = trialinfo(:, lock_col);

    % ---- accuracy + RT inclusion (identical logic to s02_epoch.m) ----
    rt_loc  = (trialinfo(:, cfg.ti.col_locresp)  - trialinfo(:, cfg.ti.col_locprobe)) / srate;
    rt_stat = (trialinfo(:, cfg.ti.col_statresp) - trialinfo(:, cfg.ti.col_statq))    / srate;

    mad_loc  = median(abs(rt_loc(~isnan(rt_loc))   - median(rt_loc,'omitnan')));
    mad_stat = median(abs(rt_stat(~isnan(rt_stat)) - median(rt_stat,'omitnan')));
    loc_degenerate  = mad_loc  < cfg.rt.degenerate_mad_thresh_s;
    stat_degenerate = mad_stat < cfg.rt.degenerate_mad_thresh_s;

    rt_ok = true(nTrials, 1);
    if cfg.rt.use_rt_filter
        if ~loc_degenerate,  rt_ok = rt_ok & local_rt_ok(rt_loc, cfg);  end
        if ~stat_degenerate, rt_ok = rt_ok & local_rt_ok(rt_stat, cfg); end
    end

    acc_ok = true(nTrials, 1);
    if cfg.correct_only
        acc_ok = (trialinfo(:, cfg.accuracy_column) == cfg.accuracy_correct_value);
    end

    boundary_samples = [];
    if isfield(EEG, 'event') && ~isempty(EEG.event)
        isBound = strcmp({EEG.event.type}, 'boundary');
        boundary_samples = [EEG.event(isBound).latency];
    end

    a = lock + win(1);
    z = lock + win(2);
    in_bounds = (a >= 1) & (z <= EEG.pnts);

    crosses_boundary = false(nTrials, 1);
    if ~isempty(boundary_samples)
        for t = 1:nTrials
            if in_bounds(t)
                crosses_boundary(t) = any(boundary_samples > a(t) & boundary_samples < z(t));
            end
        end
    end

    good = acc_ok & rt_ok & in_bounds & ~crosses_boundary;
    fprintf('  %d/%d trials good (cue-locked window) | acc_ok=%d rt_ok=%d in_bounds=%d no_boundary=%d\n', ...
        sum(good), nTrials, sum(acc_ok), sum(rt_ok), sum(in_bounds), sum(~crosses_boundary));

    if sum(good) < 10
        warning('%s: only %d good trials for alpha validation, results will be unstable. Continuing anyway.', key, sum(good));
    end

    trl = [a(good), z(good), repmat(win(1), sum(good), 1)];
    ti_good = trialinfo(good, :);

    ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
    cfg_ep = []; cfg_ep.trl = trl;
    epoched = ft_redefinetrial(cfg_ep, ft_raw);
    epoched.trialinfo = ti_good;
    elec = epoched.elec;   % capture BEFORE freqanalysis -- ft_freqdescriptives/ft_selectdata
                            % downstream do not preserve this field, so it must be
                            % stored separately rather than read back off the TFR struct.

    cfg_tf            = [];
    cfg_tf.method     = 'wavelet';
    cfg_tf.output     = 'pow';
    cfg_tf.foi        = cfg.alpha.foi;
    cfg_tf.width      = cfg.alpha.cycles;
    cfg_tf.toi        = toi;
    cfg_tf.pad        = 'nextpow2';
    cfg_tf.keeptrials = 'yes';
    tfr_all = ft_freqanalysis(cfg_tf, epoched);

    % ---- dB baseline correction: single fixed window, exact per paper
    %      (no longer an approximation -- t=0 IS cue onset for every trial) ----
    cfg_bl              = [];
    cfg_bl.baseline     = bl;
    cfg_bl.baselinetype = 'db';
    tfr_all = ft_freqbaseline(cfg_bl, tfr_all);

    ti2 = tfr_all.trialinfo;
    idx_valid_stick    = ti2(:, cfg.ti.col_val)==cfg.ti.val_valid   & ti2(:, cfg.ti.col_cue)==cfg.ti.cue_stick;
    idx_invalid_stick  = ti2(:, cfg.ti.col_val)==cfg.ti.val_invalid & ti2(:, cfg.ti.col_cue)==cfg.ti.cue_stick;
    idx_valid_avatar   = ti2(:, cfg.ti.col_val)==cfg.ti.val_valid   & ti2(:, cfg.ti.col_cue)==cfg.ti.cue_avatar;
    idx_invalid_avatar = ti2(:, cfg.ti.col_val)==cfg.ti.val_invalid & ti2(:, cfg.ti.col_cue)==cfg.ti.cue_avatar;

    n_cells = [sum(idx_valid_stick), sum(idx_invalid_stick), sum(idx_valid_avatar), sum(idx_invalid_avatar)];
    fprintf('  cell n: valid-stick=%d invalid-stick=%d valid-avatar=%d invalid-avatar=%d\n', n_cells);
    if any(n_cells < 5)
        warning('%s: at least one condition cell has <5 trials after cue-locked filtering.', key);
    end

    subj_result = struct();
    subj_result.subject            = subj;
    subj_result.key                = key;
    subj_result.n_cells            = n_cells;
    subj_result.elec               = elec;
    subj_result.tfr_valid_stick    = select_and_avg_trials(tfr_all, idx_valid_stick);
    subj_result.tfr_invalid_stick  = select_and_avg_trials(tfr_all, idx_invalid_stick);
    subj_result.tfr_valid_avatar   = select_and_avg_trials(tfr_all, idx_valid_avatar);
    subj_result.tfr_invalid_avatar = select_and_avg_trials(tfr_all, idx_invalid_avatar);

    save(out_file, 'subj_result');
    done_keys{end+1} = key; %#ok<AGROW>
    save(alpha_manifest_file, 'done_keys');
end

fprintf('\ns08_alpha_epoch_tfr complete.\n');

% ======================= local functions =======================
function ok = local_rt_ok(rt, cfg)
    ok = ~isnan(rt) & (rt >= cfg.rt.min_rt_s);
    valid_rt = rt(ok);
    if numel(valid_rt) >= 5
        med = median(valid_rt);
        madv = median(abs(valid_rt - med));
        if madv > 0
            lo = med - cfg.rt.outlier_n * madv * 1.4826;
            hi = med + cfg.rt.outlier_n * madv * 1.4826;
            ok = ok & (rt >= lo) & (rt <= hi);
        end
    end
end

function tfr_avg = select_and_avg_trials(tfr_all, idx)
    cfg_sel        = [];
    cfg_sel.trials = find(idx);
    tfr_sub = ft_selectdata(cfg_sel, tfr_all);
    cfg_avg = [];
    tfr_avg = ft_freqdescriptives(cfg_avg, tfr_sub);
end
