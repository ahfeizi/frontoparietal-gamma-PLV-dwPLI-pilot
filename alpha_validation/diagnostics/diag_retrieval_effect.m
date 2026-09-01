% s10_diagnose_retrieval_effect.m  --  DIAGNOSTIC ONLY (not part of the
% main analysis pipeline). Investigates the two significant negative
% clusters found in the RETRIEVAL window (5-6s), which the paper reports
% as NULL (all p>0.1). A significant effect exactly where the paper found
% none, combined with FAILING to replicate the paper's actual significant
% effect (encoding window), is a classic confound signature -- this
% script tries to distinguish "genuine unexpected effect" from
% "artifact/confound" rather than assuming either.
%
% WHAT THIS DOES
%   1. For each cue type's significant retrieval cluster (stick: neg
%      cluster p=0.044; avatar: neg cluster p=0.0036), extracts the
%      cluster's channel/freq/time extent from stat.negclusterslabelmat
%      and reports it in plain terms (which channels, which time window).
%   2. Plots a simple topography of the within-cluster t-statistic
%      (custom 2D scatter using elec.chanpos -- no dependency on a
%      FieldTrip layout file) so you can see by eye whether the pattern
%      looks like a plausible neural topography (fronto-parietal /
%      posterior) or a scattered/edge-channel pattern more typical of an
%      artifact.
%   3. Plots the group-average time-course of the valid-invalid alpha
%      power difference, averaged over the cluster's channels+freqs,
%      across the analysis window (not just the cluster's own time
%      extent) so you can see whether the effect is confined to the
%      pre-registered window or bleeds outward (bleeding outward, e.g.
%      ramping up toward the response, is more consistent with a
%      motor-preparation/RT confound than a discrete neural event).
%   4. For each subject, computes (a) the RT-validity effect (RT_invalid
%      - RT_valid, from s07's behavioral_by_cell.csv) and (b) the
%      cluster-averaged alpha power difference (valid-invalid, dB), then
%      tests their correlation with a permutation test (matching this
%      pipeline's existing sign-flip permutation philosophy from
%      s07/s04, rather than pulling in a parametric-only toolbox
%      function). A significant positive relationship (bigger RT effect
%      -> bigger alpha effect) is evidence FOR the RT/motor-preparation
%      confound hypothesis; no relationship does not rule it out but
%      removes this specific piece of evidence for it.
%
% OUTPUT
%   cfg.deriv_dir/alpha_diagnostics/
%     retrieval_<cue>_topo.png
%     retrieval_<cue>_timecourse.png
%     retrieval_<cue>_RTcorrelation.png
%     diagnostics_report.txt

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;

out_dir = fullfile(cfg.deriv_dir, 'alpha_diagnostics');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

S = load(fullfile(cfg.deriv_dir, 'alpha_group_stats.mat'));
all_tfr = S.all_tfr;
elec = all_tfr(1).elec;

% ---- behavioral RT-validity effect, per subject per cue type (from s07) ----
B = readtable(fullfile(cfg.deriv_dir, 'behavioral_by_cell.csv'));

rng(1);
n_perm = cfg.alpha.n_perm;

report_lines = {};
report_lines{end+1} = sprintf('Retrieval-window cluster diagnostics: %s', datestr(now));

cue_cases = struct( ...
    'label',      {'stick',                 'avatar'}, ...
    'stat',       {S.stat_retrieval_stick,   S.stat_retrieval_avatar}, ...
    'sign',       {'neg',                    'neg'}, ...
    'valid_fld',  {'tfr_valid_stick',        'tfr_valid_avatar'}, ...
    'invalid_fld',{'tfr_invalid_stick',      'tfr_invalid_avatar'});

for ci = 1:numel(cue_cases)
    cc = cue_cases(ci);
    fprintf('\n=== %s: retrieval-window %s-cluster diagnostics ===\n', cc.label, cc.sign);
    stat = cc.stat;

    clusterfield = [cc.sign 'clusters'];
    labelmatfield = [cc.sign 'clusterslabelmat'];
    if ~isfield(stat, clusterfield) || isempty(stat.(clusterfield))
        fprintf('  No %s clusters found for %s -- nothing to diagnose.\n', cc.sign, cc.label);
        report_lines{end+1} = sprintf('%s: no %s clusters found.', cc.label, cc.sign); %#ok<AGROW>
        continue
    end

    mask = stat.(labelmatfield) == 1;  % cluster #1 = smallest p, per FieldTrip convention
    p_val = stat.(clusterfield)(1).prob;

    % ---- extent ----
    chan_idx = any(any(mask,2),3);
    freq_idx = squeeze(any(any(mask,1),3));
    time_idx = squeeze(any(any(mask,1),2));

    cluster_chans = stat.label(chan_idx);
    freq_range = stat.freq(freq_idx);
    time_range = stat.time(time_idx);

    fprintf('  p = %.4f\n', p_val);
    fprintf('  channels (%d): %s\n', numel(cluster_chans), strjoin(cluster_chans, ', '));
    fprintf('  freq extent: %.1f-%.1f Hz\n', min(freq_range), max(freq_range));
    fprintf('  time extent: %.2f-%.2f s\n', min(time_range), max(time_range));

    report_lines{end+1} = sprintf('--- %s (%s cluster, p=%.4f) ---', cc.label, cc.sign, p_val); %#ok<AGROW>
    report_lines{end+1} = sprintf('  channels (%d): %s', numel(cluster_chans), strjoin(cluster_chans, ', ')); %#ok<AGROW>
    report_lines{end+1} = sprintf('  freq extent: %.1f-%.1f Hz | time extent: %.2f-%.2f s', ...
        min(freq_range), max(freq_range), min(time_range), max(time_range)); %#ok<AGROW>

    % ---- (2) topography: mean within-cluster t-value per channel ----
    t_topo = nan(numel(stat.label), 1);
    for c = find(chan_idx)'
        vals = squeeze(stat.stat(c,:,:));
        m = squeeze(mask(c,:,:));
        t_topo(c) = mean(vals(m), 'omitnan');
    end
    plot_topo(elec, stat.label, t_topo, chan_idx, ...
        sprintf('%s retrieval %s-cluster (p=%.4f): mean t-value', cc.label, cc.sign, p_val), ...
        fullfile(out_dir, sprintf('retrieval_%s_topo.png', cc.label)));

    % ---- (3) time-course of valid-invalid diff, cluster chans+freqs, full window ----
    n_sub = numel(all_tfr);
    tfr_ref = all_tfr(1).(cc.valid_fld);
    freq_mask_full = tfr_ref.freq >= min(freq_range) & tfr_ref.freq <= max(freq_range);

    diff_ts = nan(n_sub, numel(tfr_ref.time));
    for s = 1:n_sub
        v  = all_tfr(s).(cc.valid_fld);
        iv = all_tfr(s).(cc.invalid_fld);
        % per-subject channel mapping -- a subject with a manually removed
        % bad channel (e.g. sub-48 / P1) has a different channel count/order
        % than subject 1, so indices must be resolved per-subject.
        [tf_s, chan_map_s] = ismember(cluster_chans, v.label);
        if ~any(tf_s), continue; end
        chan_map_s = chan_map_s(tf_s);
        d = v.powspctrm - iv.powspctrm; % valid - invalid, dB
        diff_ts(s,:) = squeeze(mean(mean(d(chan_map_s, freq_mask_full, :), 1, 'omitnan'), 2, 'omitnan'));
    end
    plot_timecourse(tfr_ref.time, diff_ts, [min(time_range) max(time_range)], ...
        sprintf('%s retrieval %s-cluster: valid-invalid alpha power (dB)', cc.label, cc.sign), ...
        fullfile(out_dir, sprintf('retrieval_%s_timecourse.png', cc.label)));

    % ---- (4) per-subject cluster-averaged alpha diff vs RT-validity effect ----
    time_mask_full = tfr_ref.time >= min(time_range) & tfr_ref.time <= max(time_range);
    alpha_diff_subj = nan(n_sub, 1);
    subj_ids = cell(n_sub, 1);
    for s = 1:n_sub
        v  = all_tfr(s).(cc.valid_fld);
        iv = all_tfr(s).(cc.invalid_fld);
        [tf_s, chan_map_s] = ismember(cluster_chans, v.label);
        if any(tf_s)
            chan_map_s = chan_map_s(tf_s);
            d = v.powspctrm - iv.powspctrm;
            alpha_diff_subj(s) = mean(d(chan_map_s, freq_mask_full, time_mask_full), 'all', 'omitnan');
        end
        subj_ids{s} = all_tfr(s).subject;
    end

    rt_diff_subj = nan(n_sub, 1);
    for s = 1:n_sub
        r_valid   = B(strcmp(B.subject, subj_ids{s}) & strcmp(B.cue, cc.label) & strcmp(B.valid,'valid'), :);
        r_invalid = B(strcmp(B.subject, subj_ids{s}) & strcmp(B.cue, cc.label) & strcmp(B.valid,'invalid'), :);
        if ~isempty(r_valid) && ~isempty(r_invalid)
            rt_diff_subj(s) = r_invalid.rt_loc_s(1) - r_valid.rt_loc_s(1); % invalid - valid, seconds
        end
    end

    ok = ~isnan(alpha_diff_subj) & ~isnan(rt_diff_subj);
    [r_obs, p_perm] = perm_corr_test(rt_diff_subj(ok), alpha_diff_subj(ok), n_perm);
    fprintf('  RT-validity effect vs alpha diff: r=%.3f, permutation p=%.4f (N=%d, n_perm=%d)\n', ...
        r_obs, p_perm, sum(ok), n_perm);
    report_lines{end+1} = sprintf('  RT-alpha correlation: r=%.3f, perm p=%.4f (N=%d)', r_obs, p_perm, sum(ok)); %#ok<AGROW>

    plot_rt_correlation(rt_diff_subj(ok), alpha_diff_subj(ok), r_obs, p_perm, cc.label, ...
        fullfile(out_dir, sprintf('retrieval_%s_RTcorrelation.png', cc.label)));
end

report_file = fullfile(out_dir, 'diagnostics_report.txt');
fid = fopen(report_file, 'w');
fprintf(fid, '%s\n', report_lines{:});
fclose(fid);
fprintf('\nAll diagnostic figures + report saved to %s\n', out_dir);

% ======================= local functions =======================

function plot_topo(elec, labels, vals, highlight_idx, ttl, outfile)
    if isfield(elec, 'chanpos')
        pos_field = elec.chanpos;
    elseif isfield(elec, 'elecpos')
        pos_field = elec.elecpos;
    elseif isfield(elec, 'pnt')
        pos_field = elec.pnt;  % older FieldTrip / eeglab2fieldtrip naming
    else
        error('elec struct has none of chanpos/elecpos/pnt -- cannot plot topography. Fields present: %s', ...
            strjoin(fieldnames(elec), ', '));
    end
    pos = pos_field(:,1:2);
    % match elec channel order to stat.label order
    [tf, loc] = ismember(labels, elec.label);
    pos_ordered = nan(numel(labels), 2);
    pos_ordered(tf,:) = pos(loc(tf),:);

    f = figure('Visible','off','Color','w','Position',[100 100 500 450]);
    valid_v = ~isnan(vals);
    scatter(pos_ordered(~highlight_idx & valid_v,1), pos_ordered(~highlight_idx & valid_v,2), ...
        40, [0.75 0.75 0.75], 'filled'); hold on;
    scatter(pos_ordered(highlight_idx,1), pos_ordered(highlight_idx,2), ...
        140, vals(highlight_idx), 'filled', 'MarkerEdgeColor','k','LineWidth',1.2);
    colormap(f, 'jet'); cb = colorbar; cb.Label.String = 't-value';
    for i = find(highlight_idx)'
        text(pos_ordered(i,1), pos_ordered(i,2), labels{i}, ...
            'FontSize', 7, 'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
    axis equal off; title(ttl, 'Interpreter','none','FontSize',10);
    saveas(f, outfile); close(f);
end

function plot_timecourse(time_vec, diff_ts, cluster_window, ttl, outfile)
    m = mean(diff_ts, 1, 'omitnan');
    sem = std(diff_ts, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(diff_ts),1));

    f = figure('Visible','off','Color','w','Position',[100 100 700 400]);
    patch([cluster_window(1) cluster_window(2) cluster_window(2) cluster_window(1)], ...
          [min(m-sem) min(m-sem) max(m+sem) max(m+sem)], ...
          [1 0.9 0.9], 'EdgeColor','none'); hold on;
    patch([time_vec, fliplr(time_vec)], [m+sem, fliplr(m-sem)], ...
          [0.7 0.7 0.9], 'EdgeColor','none', 'FaceAlpha',0.6);
    plot(time_vec, m, 'b-', 'LineWidth', 1.8);
    yline(0, 'k:');
    xlabel('Time (s, relative to trial Start)'); ylabel('Valid - Invalid alpha power (dB)');
    title(ttl, 'Interpreter','none','FontSize',10);
    xlim([time_vec(1) time_vec(end)]);
    saveas(f, outfile); close(f);
end

function plot_rt_correlation(rt_diff, alpha_diff, r_obs, p_perm, cue_label, outfile)
    f = figure('Visible','off','Color','w','Position',[100 100 500 450]);
    scatter(rt_diff*1000, alpha_diff, 50, 'filled'); hold on;
    pfit = polyfit(rt_diff, alpha_diff, 1);
    xx = linspace(min(rt_diff), max(rt_diff), 50);
    plot(xx*1000, polyval(pfit, xx), 'r-', 'LineWidth', 1.5);
    xlabel('RT_{invalid} - RT_{valid} (ms)');
    ylabel('Cluster alpha power: valid - invalid (dB)');
    title(sprintf('%s: r=%.3f, perm p=%.4f', cue_label, r_obs, p_perm), 'FontSize', 10);
    saveas(f, outfile); close(f);
end

function [r_obs, p_perm] = perm_corr_test(x, y, n_perm)
    % Pearson correlation with a permutation-based p-value (shuffle y),
    % consistent with the sign-flip permutation approach already used in
    % s04/s07 rather than depending on the Statistics Toolbox for p.
    ok = ~isnan(x) & ~isnan(y);
    x = x(ok); y = y(ok);
    r_obs = local_corr(x, y);
    n = numel(x);
    r_null = nan(n_perm, 1);
    for p = 1:n_perm
        y_perm = y(randperm(n));
        r_null(p) = local_corr(x, y_perm);
    end
    p_perm = mean(abs(r_null) >= abs(r_obs));
end

function r = local_corr(x, y)
    xc = x - mean(x); yc = y - mean(y);
    r = sum(xc.*yc) / sqrt(sum(xc.^2) * sum(yc.^2));
end
