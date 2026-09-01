% s08a_check_trial_timeline.m  --  DIAGNOSTIC (no assumptions, just measurement)
%
% PURPOSE
%   The paper reports its within-trial analysis windows in absolute trial
%   time: ~1.5s ("shift up starts"), 3.5-4s (encoding, where cue validity
%   modulates alpha/theta -- THE effect we want to replicate), 5-6s
%   (probe/retrieval). It does not say in the excerpt we have which
%   trialinfo event that t=0 corresponds to.
%
%   Rather than guess, this script computes the median (and IQR) latency,
%   in seconds, from trial Start (col_start, col 8) to every other key
%   event (Cue, Targets, Maint, LocProbe), across all subjects. Comparing
%   these empirical numbers to the paper's 1.5 / 3.5-4 / 5-6 s structure
%   should make it clear whether t=0 = Start or t=0 = Cue, and which
%   trialinfo column the "encoding period" (3.5-4s) actually corresponds
%   to (most likely Targets or Maint onset).
%
% OUTPUT
%   Console table + trial_timeline_check.csv in cfg.deriv_dir.
%   READ THE OUTPUT before setting cfg.alpha.analysis_s / stat_latency.

clear; clc;
cfg = config();

sub_dirs = dir(fullfile(cfg.bids_dir, 'sub-*'));
sub_dirs = sub_dirs([sub_dirs.isdir]);

rows = {};
for s = 1:numel(sub_dirs)
    subj = sub_dirs(s).name;
    ti_path = cfg.authors_deriv_file(subj);
    if ~exist(ti_path, 'file'), continue; end
    L = load(ti_path, 'trialinfo'); ti = L.trialinfo;

    srate = cfg.resample_hz;
    t_start   = ti(:, cfg.ti.col_start);
    t_cue     = (ti(:, cfg.ti.col_cue_sample)     - t_start) / srate;
    t_targets = (ti(:, cfg.ti.col_targets_sample) - t_start) / srate;
    t_maint   = (ti(:, cfg.ti.col_maint)          - t_start) / srate;
    t_probe   = (ti(:, cfg.ti.col_locprobe)       - t_start) / srate;

    rows(end+1,:) = { subj, ...
        median(t_cue,'omitnan'),     iqr(t_cue), ...
        median(t_targets,'omitnan'), iqr(t_targets), ...
        median(t_maint,'omitnan'),   iqr(t_maint), ...
        median(t_probe,'omitnan'),   iqr(t_probe) }; %#ok<SAGROW>
end

T = cell2table(rows, 'VariableNames', ...
    {'subject','cue_s_med','cue_s_iqr','targets_s_med','targets_s_iqr', ...
     'maint_s_med','maint_s_iqr','probe_s_med','probe_s_iqr'});
writetable(T, fullfile(cfg.deriv_dir, 'trial_timeline_check.csv'));

fprintf('\n=== Median latency (s) from trial Start, pooled across subjects ===\n');
fprintf('  Cue onset:     %.3f s  (IQR %.3f)\n', median(T.cue_s_med),     median(T.cue_s_iqr));
fprintf('  Targets onset: %.3f s  (IQR %.3f)\n', median(T.targets_s_med), median(T.targets_s_iqr));
fprintf('  Maint onset:   %.3f s  (IQR %.3f)\n', median(T.maint_s_med),   median(T.maint_s_iqr));
fprintf('  LocProbe onset:%.3f s  (IQR %.3f)\n', median(T.probe_s_med),   median(T.probe_s_iqr));

fprintf(['\nCompare the above to the paper''s reported structure:\n' ...
    '  ~1.5 s  -- "shift up starts" (cue-type comparison window: 1.5-2.5s)\n' ...
    '  3.5-4 s -- 500ms ENCODING window (validity effect -- the one we want)\n' ...
    '  5-6 s   -- probe/retrieval window\n\n' ...
    'If Cue-onset median is close to 0s and Targets/Maint median lands near\n' ...
    '3.5-4s, t=0 = trial Start and the encoding-period column is whichever\n' ...
    'of Targets/Maint falls in that range -- use THAT column''s onset as the\n' ...
    'stat_latency reference, not col_cue_sample.\n' ...
    'If instead Cue-onset median itself is near 0s AND Maint/Targets are\n' ...
    'offset from it by roughly 3.5-4s, then t=0 = Cue onset, matching our\n' ...
    'current cfg.alpha.lock_event = ''cue'' assumption directly, and\n' ...
    'stat_latency = [3.5 4.0] can be used as-is relative to cue onset.\n' ...
    'Send me the printed medians and I will set analysis_s / stat_latency\n' ...
    'precisely instead of guessing further.\n']);

function v = iqr(x)
    x = x(~isnan(x));
    if isempty(x); v = NaN; return; end
    v = prctile(x,75) - prctile(x,25);
end
