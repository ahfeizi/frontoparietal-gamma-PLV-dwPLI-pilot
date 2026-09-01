% s09_alpha_stats.m  --  Group-level alpha validation: cluster-based
% permutation test (valid vs invalid, per cue type), restricted to the
% paper's own a-priori analysis windows.
%
% v3: two separate stat runs, matching the two validity comparisons the
% paper actually reported:
%   1. PRIMARY  -- encoding window (cfg.alpha.stat_latency_encoding =
%      [3.5 4.0] s post trial-Start). Paper found a significant validity
%      effect here. This is the effect we are trying to replicate.
%   2. SECONDARY (specificity check) -- retrieval window
%      (cfg.alpha.stat_latency_retrieval = [5.0 6.0] s). Paper found NO
%      validity effect here (all p>0.1). If our pipeline finds a
%      "significant" cluster here too, that is a red flag (suggests a
%      confound/artifact rather than a genuine replication), not a bonus
%      finding.
%
% Both use FieldTrip cluster-based permutation (cfg.alpha.n_perm = 5000),
% the paper's own statistical framework -- not this pipeline's dwPLI pair
% x band FWER approach (cfg.stat.n_perm sign-flip, used in s04/s07).

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;

M = load(fullfile(cfg.deriv_dir, 'manifest.mat')); manifest = M.manifest;

all_tfr = struct('subject', {}, 'key', {}, 'n_cells', {}, 'elec', {}, ...
    'tfr_valid_stick', {}, 'tfr_invalid_stick', {}, ...
    'tfr_valid_avatar', {}, 'tfr_invalid_avatar', {});

for r = 1:height(manifest)
    key = char(manifest.key(r));
    f = fullfile(cfg.deriv_dir, [key '_alphaTFR.mat']);
    if ~exist(f, 'file')
        warning('%s: no alpha TFR file found (run s08_alpha_epoch_tfr.m first), skipping.', key);
        continue;
    end
    L = load(f);
    all_tfr(end+1) = L.subj_result; %#ok<AGROW>
end

fprintf('Loaded alpha TFR for %d subjects.\n', numel(all_tfr));
if numel(all_tfr) < 2
    error('Need at least 2 subjects with alpha TFR results to run group stats.');
end

%% ---- Descriptive check, split by window (do BEFORE trusting inferential stats) ----
fprintf('\n=== Alpha-band (%.0f-%.0f Hz) descriptives, ENCODING window [%.1f %.1f]s, N=%d ===\n', ...
    cfg.alpha.band(1), cfg.alpha.band(2), cfg.alpha.stat_latency_encoding(1), cfg.alpha.stat_latency_encoding(2), numel(all_tfr));
report_alpha_descriptives(all_tfr, cfg.alpha.band, cfg.alpha.stat_latency_encoding);

fprintf('\n=== Alpha-band (%.0f-%.0f Hz) descriptives, RETRIEVAL window [%.1f %.1f]s, N=%d ===\n', ...
    cfg.alpha.band(1), cfg.alpha.band(2), cfg.alpha.stat_latency_retrieval(1), cfg.alpha.stat_latency_retrieval(2), numel(all_tfr));
report_alpha_descriptives(all_tfr, cfg.alpha.band, cfg.alpha.stat_latency_retrieval);

%% ---- Cluster-based permutation, valid vs invalid, per cue type ----
elec = all_tfr(1).elec;   % saved separately in s08 (not present on the TFR structs themselves)
cfg_neighb          = [];
cfg_neighb.method   = 'distance';
cfg_neighb.elec     = elec;
neighbours = ft_prepare_neighbours(cfg_neighb);

fprintf('\n--- PRIMARY: encoding-window validity effect (paper: significant) ---\n');
stat_encoding_stick  = run_cluster_stat(all_tfr, 'tfr_valid_stick',  'tfr_invalid_stick',  neighbours, cfg.alpha, cfg.alpha.stat_latency_encoding);
stat_encoding_avatar = run_cluster_stat(all_tfr, 'tfr_valid_avatar', 'tfr_invalid_avatar', neighbours, cfg.alpha, cfg.alpha.stat_latency_encoding);

fprintf('\n--- SECONDARY: retrieval-window validity effect (paper: null, specificity check) ---\n');
stat_retrieval_stick  = run_cluster_stat(all_tfr, 'tfr_valid_stick',  'tfr_invalid_stick',  neighbours, cfg.alpha, cfg.alpha.stat_latency_retrieval);
stat_retrieval_avatar = run_cluster_stat(all_tfr, 'tfr_valid_avatar', 'tfr_invalid_avatar', neighbours, cfg.alpha, cfg.alpha.stat_latency_retrieval);

save(fullfile(cfg.deriv_dir, 'alpha_group_stats.mat'), ...
    'stat_encoding_stick', 'stat_encoding_avatar', ...
    'stat_retrieval_stick', 'stat_retrieval_avatar', 'all_tfr');

%% ---- Report ----
report_lines = {};
report_lines{end+1} = sprintf('Alpha validation group stats: %s', datestr(now));
report_lines{end+1} = sprintf('Lock event: %s | N subjects: %d | n_perm: %d', ...
    cfg.alpha.lock_event, numel(all_tfr), cfg.alpha.n_perm);
report_lines{end+1} = '--- Encoding window (primary; paper: significant) ---';
report_lines{end+1} = summarize_cluster_stat(stat_encoding_stick, 'stick');
report_lines{end+1} = summarize_cluster_stat(stat_encoding_avatar, 'avatar');
report_lines{end+1} = '--- Retrieval window (secondary; paper: null, all p>0.1) ---';
report_lines{end+1} = summarize_cluster_stat(stat_retrieval_stick, 'stick');
report_lines{end+1} = summarize_cluster_stat(stat_retrieval_avatar, 'avatar');

report_file = fullfile(cfg.deriv_dir, 'alpha_validation_report.txt');
fid = fopen(report_file, 'w');
fprintf(fid, '%s\n', report_lines{:});
fclose(fid);

fprintf('\n%s\n', report_lines{:});
fprintf('Report written to %s\n', report_file);

% ======================= local functions =======================
function report_alpha_descriptives(all_tfr, band, latency)
    conds = {'tfr_valid_stick','tfr_invalid_stick','tfr_valid_avatar','tfr_invalid_avatar'};
    for c = 1:numel(conds)
        vals = nan(numel(all_tfr),1);
        for s = 1:numel(all_tfr)
            tfr = all_tfr(s).(conds{c});
            fidx = tfr.freq >= band(1) & tfr.freq <= band(2);
            tidx = tfr.time >= latency(1) & tfr.time <= latency(2);
            vals(s) = mean(tfr.powspctrm(:, fidx, tidx), 'all', 'omitnan');
        end
        fprintf('  %-20s mean dB = %+.3f  (N=%d)\n', conds{c}, mean(vals,'omitnan'), sum(~isnan(vals)));
    end
end

function stat = run_cluster_stat(all_tfr, cond_a, cond_b, neighbours, alpha_cfg, latency)
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
    cfg_stat.alpha             = 0.025;
    cfg_stat.numrandomization  = alpha_cfg.n_perm;
    cfg_stat.neighbours        = neighbours;
    cfg_stat.frequency         = alpha_cfg.band;
    cfg_stat.latency           = latency;   % restrict to the paper's a-priori window
    cfg_stat.avgovertime       = 'no';      % keep time dimension in the cluster search within this window
    if isfield(alpha_cfg, 'avgoverfreq') && alpha_cfg.avgoverfreq
        cfg_stat.avgoverfreq   = 'yes';     % paper: "averaging across the frequency band" before cluster stats
                                              % (chan x time clustering, not chan x freq x time)
    end
    cfg_stat.design            = [ones(1,n), 2*ones(1,n); 1:n, 1:n];
    cfg_stat.uvar              = 2;
    cfg_stat.ivar              = 1;

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
