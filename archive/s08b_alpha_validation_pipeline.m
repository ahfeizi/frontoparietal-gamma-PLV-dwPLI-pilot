%% s08b_alpha_validation_pipeline.m
%
% GOAL
%   Reproduce the alpha-power effect (Fig./Table item "1" referenced in
%   the paper) with our own pipeline, as a validation step BEFORE trusting
%   our gamma-band dwPLI connectivity results on the same dataset.
%
% This script implements steps (a), (b), (c) exactly as specified:
%
%   (a) New, SEPARATE epoching (not the existing maintenance-locked
%       epochs used for the dwPLI pipeline):
%         - Locked to cfg.lock_event ('cue' or 'target'), decided by
%           s08a_explore_cue_markers.m.
%         - If locked to target onset, the window [0, 0.5] s post-onset
%           is treated as the encoding-window analogue.
%         - Baseline: paper uses -500 to -100 ms pre-cue. If cue onset is
%           unavailable, the nearest defensible equivalent (-500 to -100
%           ms pre-target) is used, and this deviation is logged to
%           cfg.deviation_log and printed at the end of the run.
%
%   (b) Time-frequency parameters matched to the paper for comparability:
%         - Morlet wavelet, 2-30 Hz in 1 Hz steps
%         - 3 cycles per time window (fixed width, not scaled with freq)
%         - 50 ms time steps
%         - Alpha band = 8-12 Hz (averaged)
%         - dB baseline correction
%
%   (c) Contrast: valid vs invalid, separately for each cue type
%       (avatar / stick), using CLUSTER-BASED PERMUTATION in FieldTrip
%       (5000 permutations) - i.e. the paper's own statistical framework,
%       NOT our dwPLI pipeline's pair x band FWER approach. This keeps
%       the validation apples-to-apples against the published result.
%
% DEPENDENCIES
%   FieldTrip on path. EEGLAB *_ica.mat structures converted to FieldTrip
%   raw format via eeglab2fieldtrip (ft_hastoolbox('eeglab') style helper
%   already used elsewhere in the Pilot-2 pipeline).
%
% INPUT PER SUBJECT
%   sub-XX/ProcessedData/data_ica.mat
%     - EEG (EEGLAB struct, ICA-cleaned, picard, 1-100 Hz, re-referenced,
%       POz standing in for absent CPz, M1/M2/EOG/BIP1-24 already excluded)
%     - trialinfo (17 columns); col 4 = validity (1 valid / 2 invalid),
%       col 6 = cue type (1 stick / 2 avatar), cols 16/17 = accuracy
%
% OUTPUT
%   /derivatives/alpha_validation/sub-XX_TFR.mat   (single-subject TFR)
%   /derivatives/alpha_validation/group_stats.mat  (cluster stats, both cue types)
%   /derivatives/alpha_validation/alpha_validation_report.txt (summary)

clear; clc;

%% ---- CONFIG ----
cfg_run              = [];
cfg_run.raw_root     = '/path/to/ds003702';                      % EDIT ME
cfg_run.deriv_root   = '/path/to/derivatives/alpha_validation';   % EDIT ME
cfg_run.lock_event   = 'target';   % 'cue' | 'target' -- SET AFTER s08a REPORT
cfg_run.manual_bad_channels = get_pilot2_manual_bad_channels();   % reuse existing map, e.g. sub-48 P1
cfg_run.overwrite    = false;      % resume-capable, like the rest of the Pilot-2 pipeline

if ~exist(cfg_run.deriv_root, 'dir'); mkdir(cfg_run.deriv_root); end

% --- Deviation log: filled in only if lock_event == 'target' ---
deviation_log = {};
if strcmp(cfg_run.lock_event, 'target')
    deviation_log{end+1} = sprintf([ ...
        'DEVIATION: cue-onset marker could not be reliably identified in ' ...
        'the raw .vmrk stream (see s08a_explore_cue_markers.m output). ' ...
        'Epochs are locked to target onset instead. The [0, 0.5] s ' ...
        'post-target window is used as an encoding-window analogue for ' ...
        'the paper''s cue-locked alpha window. Baseline is -500 to -100 ' ...
        'ms pre-TARGET (not pre-cue). This should be reported explicitly ' ...
        'as a methodological deviation, not silently treated as equivalent.']);
end

subjects = get_valid_subject_list(cfg_run.raw_root); % excludes sub-08/42/47

%% ---- PER-SUBJECT: EPOCHING + TFR ----
all_tfr = struct('subject', {}, 'tfr_valid_stick', {}, 'tfr_invalid_stick', {}, ...
                  'tfr_valid_avatar', {}, 'tfr_invalid_avatar', {});

manifest_file = fullfile(cfg_run.deriv_root, 'manifest.mat');
done_subjects = {};
if exist(manifest_file, 'file')
    m = load(manifest_file);
    done_subjects = m.done_subjects;
end

for s = 1:numel(subjects)
    subj = subjects{s};
    out_file = fullfile(cfg_run.deriv_root, sprintf('%s_TFR.mat', subj));

    if ismember(subj, done_subjects) && exist(out_file, 'file') && ~cfg_run.overwrite
        fprintf('%s already processed, skipping (resume).\n', subj);
        r = load(out_file);
        all_tfr(end+1) = r.subj_result; %#ok<AGROW>
        continue
    end

    fprintf('\n=== Processing %s (lock_event = %s) ===\n', subj, cfg_run.lock_event);

    ica_file = fullfile(cfg_run.raw_root, subj, 'ProcessedData', 'data_ica.mat');
    if ~exist(ica_file, 'file')
        warning('%s: data_ica.mat not found, skipping.', subj);
        continue
    end
    L = load(ica_file); % expects L.EEG (EEGLAB struct) and L.trialinfo (17 cols)
    EEG       = L.EEG;
    trialinfo = L.trialinfo;

    % Apply subject-specific manual bad-channel overrides (same mechanism
    % as the dwPLI pipeline, e.g. sub-48 parietal P1 artifact).
    if isKey(cfg_run.manual_bad_channels, subj)
        bad_ch = cfg_run.manual_bad_channels(subj);
        EEG = pop_select(EEG, 'nochannel', bad_ch);
    end

    % --- Convert to FieldTrip raw ---
    ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
    ft_raw.trialinfo = trialinfo;

    % --- (a) NEW epoch definition, independent of the maintenance-locked
    %     epochs used elsewhere in the Pilot-2 pipeline ---
    cfg_ep = [];
    switch cfg_run.lock_event
        case 'cue'
            cfg_ep.trl = build_trial_definition(EEG, trialinfo, 'cue');
            baseline_win = [-0.5 -0.1];   % pre-cue, matches paper exactly
            analysis_win = [-0.5 1.0];    % generous window covering baseline + response
        case 'target'
            cfg_ep.trl = build_trial_definition(EEG, trialinfo, 'target');
            baseline_win = [-0.5 -0.1];   % pre-target (deviation, logged above)
            analysis_win = [-0.5 0.5];    % encoding-window analogue, 0-500ms post-target
        otherwise
            error('cfg_run.lock_event must be ''cue'' or ''target''.');
    end
    epoched = ft_redefinetrial(cfg_ep, ft_raw);

    % --- (b) TF decomposition matched to the paper ---
    cfg_tf            = [];
    cfg_tf.method     = 'wavelet';
    cfg_tf.output     = 'pow';
    cfg_tf.foi        = 2:1:30;        % 2-30 Hz, 1 Hz steps
    cfg_tf.width      = 3;             % 3 cycles per time window, fixed
    cfg_tf.toi        = analysis_win(1):0.05:analysis_win(2); % 50 ms steps
    cfg_tf.pad        = 'nextpow2';
    cfg_tf.keeptrials = 'yes';
    tfr_all = ft_freqanalysis(cfg_tf, epoched);

    % --- dB baseline correction ---
    cfg_bl           = [];
    cfg_bl.baseline  = baseline_win;
    cfg_bl.baselinetype = 'db';
    tfr_all = ft_freqbaseline(cfg_bl, tfr_all);

    % --- Split by condition: validity (col 4) x cue type (col 6) ---
    ti = tfr_all.trialinfo;
    idx_valid_stick    = ti(:,4)==1 & ti(:,6)==1;
    idx_invalid_stick  = ti(:,4)==2 & ti(:,6)==1;
    idx_valid_avatar   = ti(:,4)==1 & ti(:,6)==2;
    idx_invalid_avatar = ti(:,4)==2 & ti(:,6)==2;

    subj_result = struct();
    subj_result.subject           = subj;
    subj_result.tfr_valid_stick    = select_and_avg_trials(tfr_all, idx_valid_stick);
    subj_result.tfr_invalid_stick  = select_and_avg_trials(tfr_all, idx_invalid_stick);
    subj_result.tfr_valid_avatar   = select_and_avg_trials(tfr_all, idx_valid_avatar);
    subj_result.tfr_invalid_avatar = select_and_avg_trials(tfr_all, idx_invalid_avatar);

    save(out_file, 'subj_result');
    all_tfr(end+1) = subj_result; %#ok<AGROW>

    done_subjects{end+1} = subj; %#ok<AGROW>
    save(manifest_file, 'done_subjects'); % incremental write, survives crashes
end

%% ---- ALPHA-BAND SANITY CHECK (descriptive, pre-registration style) ----
% Quick check that we replicate the qualitative alpha pattern before
% running inferential stats: average 8-12 Hz power in the analysis
% window, per subject, per condition.
fprintf('\n=== Alpha-band (8-12 Hz) descriptive summary ===\n');
report_alpha_descriptives(all_tfr);

%% ---- (c) GROUP-LEVEL CLUSTER-BASED PERMUTATION (valid vs invalid) ----
% Run separately per cue type, using the PAPER'S statistical framework
% (cluster-based permutation, FieldTrip, 5000 permutations) rather than
% our own dwPLI pipeline's pair x band FWER correction, so the
% comparison to the published effect is apples-to-apples.

elec = all_tfr(1).tfr_valid_stick.elec; % assumes consistent electrode layout across subjects
cfg_neighb          = [];
cfg_neighb.method   = 'distance';
cfg_neighb.elec     = elec;
neighbours = ft_prepare_neighbours(cfg_neighb);

stat_stick  = run_cluster_stat(all_tfr, 'tfr_valid_stick',  'tfr_invalid_stick',  neighbours);
stat_avatar = run_cluster_stat(all_tfr, 'tfr_valid_avatar', 'tfr_invalid_avatar', neighbours);

save(fullfile(cfg_run.deriv_root, 'group_stats.mat'), 'stat_stick', 'stat_avatar');

%% ---- REPORT ----
report_lines = {};
report_lines{end+1} = sprintf('Alpha validation run: %s', datestr(now));
report_lines{end+1} = sprintf('Lock event: %s', cfg_run.lock_event);
report_lines{end+1} = sprintf('N subjects processed: %d', numel(all_tfr));
report_lines = [report_lines, deviation_log];
report_lines{end+1} = summarize_cluster_stat(stat_stick, 'stick');
report_lines{end+1} = summarize_cluster_stat(stat_avatar, 'avatar');

report_file = fullfile(cfg_run.deriv_root, 'alpha_validation_report.txt');
fid = fopen(report_file, 'w');
fprintf(fid, '%s\n', report_lines{:});
fclose(fid);

fprintf('\nDone. Report written to %s\n', report_file);

%% ---- Local helper functions ----

function bad_map = get_pilot2_manual_bad_channels()
    bad_map = containers.Map('KeyType', 'char', 'ValueType', 'any');
    bad_map('sub-48') = {'P1'}; % block-localized artifact, std=153uV vs ~5-7uV neighbors
end

function subs = get_valid_subject_list(raw_root)
    d = dir(fullfile(raw_root, 'sub-*'));
    excluded = {'sub-08', 'sub-42', 'sub-47'};
    subs = setdiff({d([d.isdir]).name}, excluded);
end

function trl = build_trial_definition(EEG, trialinfo, lock_event)
    % Builds an FieldTrip-style trl matrix [begsample endsample offset]
    % locked either to the cue-onset marker or the target-onset marker,
    % using the trial-by-trial sample indices already present in the
    % EEGLAB EEG.event structure. The exact event 'type' string to match
    % must be confirmed against s08a_explore_cue_markers.m output before
    % running this in production; placeholders below assume the codes
    % identified as most consistent pre-maintenance markers.
    fs = EEG.srate; % 500 Hz
    switch lock_event
        case 'cue'
            evtype_stick  = 'S ??';  % FILL IN from s08a RECOMMENDATION
            evtype_avatar = 'S ??';  % FILL IN from s08a RECOMMENDATION
            prestim  = round(0.5 * fs);
            poststim = round(1.0 * fs);
        case 'target'
            evtype_stick  = 'S ??';  % FILL IN from s08a RECOMMENDATION
            evtype_avatar = 'S ??';  % FILL IN from s08a RECOMMENDATION
            prestim  = round(0.5 * fs);
            poststim = round(0.5 * fs);
    end

    onset_types = {evtype_stick, evtype_avatar};
    is_onset = ismember({EEG.event.type}, onset_types);
    onset_samples = round([EEG.event(is_onset).latency]);

    n_trials = min(numel(onset_samples), size(trialinfo, 1));
    trl = zeros(n_trials, 3);
    for t = 1:n_trials
        trl(t,1) = onset_samples(t) - prestim;
        trl(t,2) = onset_samples(t) + poststim;
        trl(t,3) = -prestim;
    end
    trl = trl(trl(:,1) > 0, :); % drop trials that would run before recording start
end

function tfr_avg = select_and_avg_trials(tfr_all, idx)
    cfg_sel        = [];
    cfg_sel.trials = find(idx);
    tfr_sub = ft_selectdata(cfg_sel, tfr_all);
    cfg_avg = [];
    tfr_avg = ft_freqdescriptives(cfg_avg, tfr_sub);
end

function report_alpha_descriptives(all_tfr)
    alpha_band = [8 12];
    conds = {'tfr_valid_stick','tfr_invalid_stick','tfr_valid_avatar','tfr_invalid_avatar'};
    for c = 1:numel(conds)
        vals = nan(numel(all_tfr),1);
        for s = 1:numel(all_tfr)
            tfr = all_tfr(s).(conds{c});
            fidx = tfr.freq >= alpha_band(1) & tfr.freq <= alpha_band(2);
            vals(s) = mean(tfr.powspctrm(:, fidx, :), 'all', 'omitnan');
        end
        fprintf('  %-20s mean(8-12Hz dB) = %+.3f  (N=%d)\n', conds{c}, mean(vals,'omitnan'), sum(~isnan(vals)));
    end
end

function stat = run_cluster_stat(all_tfr, cond_a, cond_b, neighbours)
    n = numel(all_tfr);
    tfr_a = cell(1,n); tfr_b = cell(1,n);
    for s = 1:n
        tfr_a{s} = all_tfr(s).(cond_a);
        tfr_b{s} = all_tfr(s).(cond_b);
    end

    cfg_stat                  = [];
    cfg_stat.method            = 'montecarlo';
    cfg_stat.statistic         = 'ft_statfun_depsamplesT';
    cfg_stat.correctm          = 'cluster';
    cfg_stat.clusteralpha      = 0.05;
    cfg_stat.clusterstatistic  = 'maxsum';
    cfg_stat.minnbchan         = 2;
    cfg_stat.tail              = 0;
    cfg_stat.clustertail       = 0;
    cfg_stat.alpha             = 0.025;   % two-tailed
    cfg_stat.numrandomization  = 5000;    % matches paper
    cfg_stat.neighbours        = neighbours;
    cfg_stat.frequency         = [8 12];  % alpha band, matches paper's contrast band
    cfg_stat.design            = [ones(1,n), 2*ones(1,n); 1:n, 1:n];
    cfg_stat.uvar               = 2;
    cfg_stat.ivar               = 1;

    stat = ft_freqstatistics(cfg_stat, tfr_a{:}, tfr_b{:});
end

function line = summarize_cluster_stat(stat, label)
    if isfield(stat, 'posclusters') && ~isempty(stat.posclusters)
        p_pos = stat.posclusters(1).prob;
    else
        p_pos = NaN;
    end
    if isfield(stat, 'negclusters') && ~isempty(stat.negclusters)
        p_neg = stat.negclusters(1).prob;
    else
        p_neg = NaN;
    end
    line = sprintf('Cue type %-6s: smallest pos-cluster p=%.4f, smallest neg-cluster p=%.4f', ...
        label, p_pos, p_neg);
end
