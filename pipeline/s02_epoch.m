% s02_epoch.m  --  Build the maintenance-window trial list for each subject
% (Pilot 2 / ds003702), using the AUTHORS' validated trialinfo instead of
% re-parsing raw marker pairs ourselves (see config.m note).
%
% For each subject in the manifest:
%   1. Load our own cleaned continuous .set (from s01_preprocess.m)
%   2. Load trialinfo from the authors' data_ica.mat (ONLY trialinfo,
%      not their filtered EEG data)
%   3. For each trial, define the maintenance-window lock sample from
%      trialinfo column col_maint, and compute RT for both responses
%   4. Mark a trial 'good' only if: accuracy correct (if cfg.correct_only),
%      RT passes both filters (if cfg.rt.use_rt_filter), the maintenance
%      window stays within the recording, and the window does NOT cross
%      any 'boundary' event (VR-headset-removal break / ASR discontinuity)
%   5. Save <key>_trials.mat, matching Pilot 1's s02 output convention

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

srate = cfg.resample_hz;
win = round(cfg.maint_window * srate);   % [start end] samples rel. to Maint

% Cross-subject QC log for the automated degenerate-RT check (see config.m
% cfg.rt.degenerate_mad_thresh_s). Populated as we go, written once at the end.
qc_rows = {};

for r = 1:height(manifest)
    key  = char(manifest.key(r));
    subj = char(manifest.subject(r));
    fprintf('\n=== %s ===\n', key);

    % ---- our cleaned continuous data (for pnts/boundary info only) ----
    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);

    % ---- authors' trialinfo (NOT their filtered EEG data) ----
    ti_path = cfg.authors_deriv_file(subj);
    if ~exist(ti_path, 'file')
        warning('%s: trialinfo file not found at %s, skipping.', subj, ti_path);
        continue;
    end
    TI = load(ti_path, 'trialinfo');
    trialinfo = TI.trialinfo;
    nTrials = size(trialinfo, 1);

    % ---- boundary event sample positions (VR-break / ASR discontinuities) ----
    boundary_samples = [];
    if isfield(EEG, 'event') && ~isempty(EEG.event)
        isBound = strcmp({EEG.event.type}, 'boundary');
        boundary_samples = [EEG.event(isBound).latency];
    end
    fprintf('  %d boundary events in this subject''s cleaned data\n', numel(boundary_samples));

    % ---- pull relevant columns ----
    lock    = trialinfo(:, cfg.ti.col_maint);
    cueV    = trialinfo(:, cfg.ti.col_cue);
    validV  = trialinfo(:, cfg.ti.col_val);
    locV    = trialinfo(:, cfg.ti.col_loc);
    accLoc  = trialinfo(:, cfg.ti.col_acc_loc);
    accStat = trialinfo(:, cfg.ti.col_acc_stat);

    rt_loc  = (trialinfo(:, cfg.ti.col_locresp)  - trialinfo(:, cfg.ti.col_locprobe)) / srate;
    rt_stat = (trialinfo(:, cfg.ti.col_statresp) - trialinfo(:, cfg.ti.col_statq))    / srate;

    fprintf('  RT_loc  median=%.3f s, range=[%.3f %.3f]\n', median(rt_loc,'omitnan'), min(rt_loc), max(rt_loc));
    fprintf('  RT_stat median=%.3f s, range=[%.3f %.3f]\n', median(rt_stat,'omitnan'), min(rt_stat), max(rt_stat));

    % ---- automated degenerate-RT detection (per subject, per response type) ----
    % Robust spread (MAD) below cfg.rt.degenerate_mad_thresh_s means the
    % "RT" is not a genuine self-paced response (task design has no
    % response-window cut-off -- Gregory, Wang & Kessler 2022), so filtering
    % on it would be meaningless. Flag and disable that filter for this
    % subject/response-type only; the OTHER RT type and accuracy still apply.
    mad_loc  = median(abs(rt_loc(~isnan(rt_loc))   - median(rt_loc,'omitnan')));
    mad_stat = median(abs(rt_stat(~isnan(rt_stat)) - median(rt_stat,'omitnan')));
    loc_degenerate  = mad_loc  < cfg.rt.degenerate_mad_thresh_s;
    stat_degenerate = mad_stat < cfg.rt.degenerate_mad_thresh_s;
    if loc_degenerate
        fprintf('  [QC FLAG] RT_loc MAD = %.4f s (< %.3f s threshold) -- treating as degenerate, RT_loc filter DISABLED for %s\n', ...
            mad_loc, cfg.rt.degenerate_mad_thresh_s, subj);
    end
    if stat_degenerate
        fprintf('  [QC FLAG] RT_stat MAD = %.4f s (< %.3f s threshold) -- treating as degenerate, RT_stat filter DISABLED for %s\n', ...
            mad_stat, cfg.rt.degenerate_mad_thresh_s, subj);
    end
    qc_rows(end+1,:) = {subj, mad_loc, mad_stat, loc_degenerate, stat_degenerate}; %#ok<SAGROW>

    % ---- RT outlier bounds (per-subject, MAD-based; skip degenerate types) ----
    rt_ok = true(nTrials, 1);
    if cfg.rt.use_rt_filter
        if ~loc_degenerate,  rt_ok = rt_ok & local_rt_ok(rt_loc, cfg);  end
        if ~stat_degenerate, rt_ok = rt_ok & local_rt_ok(rt_stat, cfg); end
    end

    % ---- assemble trials struct ----
    trials = struct('lock_sample',{}, 'cue',{}, 'valid',{}, 'loc',{}, ...
        'acc_loc',{}, 'acc_stat',{}, 'rt_loc',{}, 'rt_stat',{}, 'good',{}, 'reason',{});

    n_excl_acc = 0; n_excl_rt = 0; n_excl_bounds = 0; n_excl_boundary = 0;

    for t = 1:nTrials
        ls = lock(t);
        a = ls + win(1); z = ls + win(2);
        reason = '';

        in_bounds = (a >= 1) && (z <= EEG.pnts);
        if ~in_bounds
            reason = 'out_of_bounds'; n_excl_bounds = n_excl_bounds + 1;
        end

        crosses_boundary = false;
        if in_bounds && ~isempty(boundary_samples)
            crosses_boundary = any(boundary_samples > a & boundary_samples < z);
            if crosses_boundary
                reason = 'crosses_boundary'; n_excl_boundary = n_excl_boundary + 1;
            end
        end

        acc_ok = true;
        if cfg.correct_only
            acc_ok = (trialinfo(t, cfg.accuracy_column) == cfg.accuracy_correct_value);
            if ~acc_ok && isempty(reason)
                reason = 'incorrect'; n_excl_acc = n_excl_acc + 1;
            end
        end

        rt_this_ok = rt_ok(t);
        if ~rt_this_ok && isempty(reason)
            reason = 'rt_outlier_or_anticipatory'; n_excl_rt = n_excl_rt + 1;
        end

        good = in_bounds && ~crosses_boundary && acc_ok && rt_this_ok;

        trials(end+1) = struct('lock_sample', ls, 'cue', cueV(t), 'valid', validV(t), ...
            'loc', locV(t), 'acc_loc', accLoc(t), 'acc_stat', accStat(t), ...
            'rt_loc', rt_loc(t), 'rt_stat', rt_stat(t), 'good', good, 'reason', reason); %#ok<SAGROW>
    end

    nGood = sum([trials.good]);
    nSocialGood    = sum([trials.good] & [trials.cue]==cfg.ti.cue_avatar);
    nNonsocialGood = sum([trials.good] & [trials.cue]==cfg.ti.cue_stick);
    nValidGood     = sum([trials.good] & [trials.valid]==cfg.ti.val_valid);
    nInvalidGood   = sum([trials.good] & [trials.valid]==cfg.ti.val_invalid);

    fprintf('  %d/%d trials good | excluded: %d out-of-bounds, %d cross-boundary, %d incorrect, %d RT\n', ...
        nGood, nTrials, n_excl_bounds, n_excl_boundary, n_excl_acc, n_excl_rt);
    fprintf('  good breakdown: avatar=%d stick=%d | valid=%d invalid=%d\n', ...
        nSocialGood, nNonsocialGood, nValidGood, nInvalidGood);

    save(fullfile(cfg.deriv_dir, [key '_trials.mat']), 'trials');
end

fprintf('\ns02_epoch complete.\n');

% ---- write cross-subject RT QC flags ----
QC = cell2table(qc_rows, 'VariableNames', ...
    {'subject','mad_rt_loc_s','mad_rt_stat_s','rt_loc_degenerate','rt_stat_degenerate'});
writetable(QC, fullfile(cfg.deriv_dir, 'rt_qc_flags.csv'));
nFlaggedLoc  = sum(QC.rt_loc_degenerate);
nFlaggedStat = sum(QC.rt_stat_degenerate);
fprintf('RT QC: %d/%d subjects flagged degenerate RT_loc, %d/%d flagged degenerate RT_stat.\n', ...
    nFlaggedLoc, height(QC), nFlaggedStat, height(QC));
fprintf('See %s for the full per-subject breakdown.\n', fullfile(cfg.deriv_dir, 'rt_qc_flags.csv'));

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
