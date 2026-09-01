% s11_diagnose_encoding_topo.m  --  DIAGNOSTIC. Compares the SPATIAL
% pattern of our encoding-window clusters against the paper's reported
% channels, not just the p-values.
%
% Paper's exact reported significant posterior alpha clusters (valid vs
% invalid, encoding window), from the Fig. 4 caption (Gregory, Wang &
% Kessler, 2022, SCAN):
%   Social (avatar) cue : CP1, P3, Pz, POz, P1, P2, PO3, Oz
%   Non-social (stick) cue: Pz, POz, O1, O2, PO5, PO3, PO4, PO6, Oz
%
% The paper also defines broader anterior/posterior ROI groupings used
% for its TFR display panels (not itself a statistical test, but useful
% here to classify whether our secondary POSITIVE cluster -- which the
% paper does not explicitly report -- looks anterior, posterior, or
% mixed):
%   Anterior : Fp1,Fpz,Fp2,AF7,AF3,AF4,AF8,F7,F5,F3,F1,Fz,F2,F4,F6,F8,
%              FT7,FC5,FC3,FC1,FCz,FC2,FC4,FC6,FT8
%   Posterior: TP7,CP5,CP3,CP1,CP2,CP4,CP6,TP8,P7,P5,P3,P1,Pz,P2,P4,P6,P8,
%              PO7,PO5,PO3,POz,PO4,PO6,PO8,O1,Oz,O2
%
% For each cue type's PRIMARY negative cluster (the replicated effect,
% p=0.0002 both cue types) and SECONDARY positive cluster (trend, not
% reported in the paper), this script:
%   1. Reports channel-list overlap with the paper's specific
%      significant-cluster channels (Jaccard-style: how many of the
%      paper's channels are in ours, and vice versa).
%   2. Classifies cluster channels as anterior/posterior/neither using
%      the paper's broader ROI split.
%   3. Plots topography with paper-reported channels marked distinctly
%      (filled diamond) vs. our-cluster-only channels (circle) for
%      visual comparison.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;

out_dir = fullfile(cfg.deriv_dir, 'alpha_diagnostics');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

S = load(fullfile(cfg.deriv_dir, 'alpha_group_stats.mat'));
all_tfr = S.all_tfr;
elec = all_tfr(1).elec;

% ---- paper's reference channel sets ----
paper_sig = struct( ...
    'avatar', {{'CP1','P3','Pz','POz','P1','P2','PO3','Oz'}}, ...
    'stick',  {{'Pz','POz','O1','O2','PO5','PO3','PO4','PO6','Oz'}});

paper_anterior = {'Fp1','Fpz','Fp2','AF7','AF3','AF4','AF8','F7','F5','F3','F1','Fz','F2','F4','F6','F8', ...
    'FT7','FC5','FC3','FC1','FCz','FC2','FC4','FC6','FT8'};
paper_posterior = {'TP7','CP5','CP3','CP1','CP2','CP4','CP6','TP8','P7','P5','P3','P1','Pz','P2','P4','P6','P8', ...
    'PO7','PO5','PO3','POz','PO4','PO6','PO8','O1','Oz','O2'};

cue_cases = struct( ...
    'label', {'stick', 'avatar'}, ...
    'stat',  {S.stat_encoding_stick, S.stat_encoding_avatar});

report_lines = {};
report_lines{end+1} = sprintf('Encoding-window cluster topography comparison: %s', datestr(now));

for ci = 1:numel(cue_cases)
    cc = cue_cases(ci);
    stat = cc.stat;
    fprintf('\n=== %s: encoding-window cluster topography vs. paper ===\n', cc.label);
    report_lines{end+1} = sprintf('--- %s ---', cc.label); %#ok<AGROW>

    cluster_time_ranges = struct('neg', [], 'pos', []);
    cluster_chans_by_sign = struct('neg', {{}}, 'pos', {{}});

    for signtype = {'neg','pos'}
        sgn = signtype{1};
        clusterfield  = [sgn 'clusters'];
        labelmatfield = [sgn 'clusterslabelmat'];

        if ~isfield(stat, clusterfield) || isempty(stat.(clusterfield))
            fprintf('  [%s] no clusters found.\n', sgn);
            continue
        end

        mask = stat.(labelmatfield) == 1;
        p_val = stat.(clusterfield)(1).prob;
        chan_idx = any(any(mask,2),3);
        cluster_chans = stat.label(chan_idx);

        % ---- time extent (dims may be chan x time if avgoverfreq='yes',
        %      or chan x freq x time otherwise -- handle both) ----
        if ndims(mask) == 2
            time_idx = any(mask, 1);
        else
            time_idx = squeeze(any(any(mask,1),2));
        end
        time_range = stat.time(time_idx);

        ref_chans = paper_sig.(cc.label);
        n_overlap = numel(intersect(cluster_chans, ref_chans));
        pct_of_paper_covered  = 100 * n_overlap / numel(ref_chans);
        pct_of_ours_in_paper  = 100 * n_overlap / numel(cluster_chans);

        n_ant = numel(intersect(cluster_chans, paper_anterior));
        n_post = numel(intersect(cluster_chans, paper_posterior));

        fprintf('  [%s, p=%.4f] %d channels, time extent %.2f-%.2fs: %s\n', ...
            sgn, p_val, numel(cluster_chans), min(time_range), max(time_range), strjoin(cluster_chans, ', '));
        fprintf('    overlap with paper''s reported cluster (%s): %d/%d paper channels covered (%.0f%%), %d/%d of ours are in paper''s list (%.0f%%)\n', ...
            strjoin(ref_chans, ','), n_overlap, numel(ref_chans), pct_of_paper_covered, n_overlap, numel(cluster_chans), pct_of_ours_in_paper);
        fprintf('    anterior/posterior split (paper ROI def.): %d anterior, %d posterior, %d neither\n', ...
            n_ant, n_post, numel(cluster_chans)-n_ant-n_post);

        report_lines{end+1} = sprintf('  [%s, p=%.4f] %d chans, time %.2f-%.2fs | paper overlap: %d/%d paper covered (%.0f%%), %d/%d ours-in-paper (%.0f%%) | ant=%d post=%d neither=%d', ...
            sgn, p_val, numel(cluster_chans), min(time_range), max(time_range), n_overlap, numel(ref_chans), pct_of_paper_covered, ...
            n_overlap, numel(cluster_chans), pct_of_ours_in_paper, n_ant, n_post, numel(cluster_chans)-n_ant-n_post); %#ok<AGROW>

        plot_topo_vs_paper(elec, stat.label, chan_idx, ref_chans, ...
            sprintf('%s encoding %s-cluster (p=%.4f) vs paper channels', cc.label, sgn, p_val), ...
            fullfile(out_dir, sprintf('encoding_%s_%s_topo_vs_paper.png', cc.label, sgn)));

        cluster_time_ranges.(sgn) = [min(time_range) max(time_range)];
        cluster_chans_by_sign.(sgn) = cluster_chans;
    end

    % ---- combined timeline: are neg and pos clusters sequential or overlapping? ----
    if ~isempty(cluster_time_ranges.neg) && ~isempty(cluster_time_ranges.pos)
        overlap_start = max(cluster_time_ranges.neg(1), cluster_time_ranges.pos(1));
        overlap_end   = min(cluster_time_ranges.neg(2), cluster_time_ranges.pos(2));
        is_overlapping = overlap_start < overlap_end;
        if is_overlapping
            fprintf('  --> neg and pos clusters OVERLAP in time (%.2f-%.2fs shared).\n', overlap_start, overlap_end);
            report_lines{end+1} = sprintf('  --> neg/pos clusters OVERLAP in time (%.2f-%.2fs shared).', overlap_start, overlap_end); %#ok<AGROW>
        else
            gap = max(cluster_time_ranges.pos(1) - cluster_time_ranges.neg(2), cluster_time_ranges.neg(1) - cluster_time_ranges.pos(2));
            fprintf('  --> neg and pos clusters are SEQUENTIAL (gap = %.2fs) -- consistent with an early desynchronization followed by a later rebound.\n', gap);
            report_lines{end+1} = sprintf('  --> neg/pos clusters SEQUENTIAL (gap = %.2fs).', gap); %#ok<AGROW>
        end

        plot_cluster_timeline(all_tfr, cc.label, cluster_chans_by_sign, cluster_time_ranges, ...
            fullfile(out_dir, sprintf('encoding_%s_cluster_timeline.png', cc.label)));
    end
end

report_file = fullfile(out_dir, 'encoding_topo_report.txt');
fid = fopen(report_file, 'w');
fprintf(fid, '%s\n', report_lines{:});
fclose(fid);
fprintf('\nReport + figures saved to %s\n', out_dir);

% ======================= local functions =======================
function plot_topo_vs_paper(elec, labels, our_idx, ref_chans, ttl, outfile)
    if isfield(elec, 'chanpos'), pos_field = elec.chanpos;
    elseif isfield(elec, 'elecpos'), pos_field = elec.elecpos;
    elseif isfield(elec, 'pnt'), pos_field = elec.pnt;
    else, error('elec struct has no recognizable position field.');
    end
    pos = pos_field(:,1:2);
    [tf, loc] = ismember(labels, elec.label);
    pos_ordered = nan(numel(labels), 2);
    pos_ordered(tf,:) = pos(loc(tf),:);

    ref_idx = ismember(labels, ref_chans);
    both_idx   = our_idx & ref_idx;
    ours_only  = our_idx & ~ref_idx;
    paper_only = ~our_idx & ref_idx;
    neither    = ~our_idx & ~ref_idx;

    f = figure('Visible','off','Color','w','Position',[100 100 550 500]);
    scatter(pos_ordered(neither,1), pos_ordered(neither,2), 30, [0.85 0.85 0.85], 'filled'); hold on;
    scatter(pos_ordered(ours_only,1), pos_ordered(ours_only,2), 100, [0.2 0.4 0.9], 'o', 'filled', 'MarkerEdgeColor','k');
    scatter(pos_ordered(paper_only,1), pos_ordered(paper_only,2), 100, [0.9 0.4 0.1], 'd', 'filled', 'MarkerEdgeColor','k');
    scatter(pos_ordered(both_idx,1), pos_ordered(both_idx,2), 130, [0.1 0.7 0.2], 'p', 'filled', 'MarkerEdgeColor','k');
    legend({'neither','our cluster only','paper only','both (match)'}, 'Location','southoutside','Orientation','horizontal','FontSize',7);
    for i = find(our_idx | ref_idx)'
        text(pos_ordered(i,1), pos_ordered(i,2), labels{i}, 'FontSize',6, 'HorizontalAlignment','center','VerticalAlignment','bottom');
    end
    axis equal off; title(ttl, 'Interpreter','none','FontSize',9);
    saveas(f, outfile); close(f);
end

function plot_cluster_timeline(all_tfr, cue_label, cluster_chans_by_sign, cluster_time_ranges, outfile)
    % Shows the group-average valid-invalid alpha power time-course
    % (averaged separately over the neg-cluster's own channels and the
    % pos-cluster's own channels), with each cluster's exact time extent
    % shaded, over a window spanning a bit before/after the encoding
    % period so the sequence (or overlap) is visible in context.
    valid_fld   = sprintf('tfr_valid_%s', cue_label);
    invalid_fld = sprintf('tfr_invalid_%s', cue_label);
    tfr_ref = all_tfr(1).(valid_fld);
    freq_mask = tfr_ref.freq >= 8 & tfr_ref.freq <= 12;   % alpha band
    win_mask  = tfr_ref.time >= 3.0 & tfr_ref.time <= 4.5; % context window around encoding (3.5-4.0s)
    t = tfr_ref.time(win_mask);

    n_sub = numel(all_tfr);
    colors = struct('neg',[0.1 0.3 0.8], 'pos',[0.8 0.3 0.1]);

    f = figure('Visible','off','Color','w','Position',[100 100 750 420]); hold on;
    y_all = [];
    for signtype = {'neg','pos'}
        sgn = signtype{1};
        chans = cluster_chans_by_sign.(sgn);
        if isempty(chans), continue; end

        diff_ts = nan(n_sub, sum(win_mask));
        for s = 1:n_sub
            v  = all_tfr(s).(valid_fld);
            iv = all_tfr(s).(invalid_fld);
            % per-subject channel mapping -- subjects with a manually
            % removed bad channel (e.g. sub-48 / P1) have a different
            % channel count/order than subject 1, so chan indices must
            % be resolved against THIS subject's own label list, not a
            % shared one computed once from subject 1.
            [tf_s, chan_map_s] = ismember(chans, v.label);
            if ~any(tf_s)
                continue  % none of this cluster's channels exist for this subject -- skip
            end
            chan_map_s = chan_map_s(tf_s);
            d = v.powspctrm - iv.powspctrm;
            diff_ts(s,:) = squeeze(mean(mean(d(chan_map_s, freq_mask, win_mask), 1, 'omitnan'), 2, 'omitnan'));
        end
        m = mean(diff_ts, 1, 'omitnan');
        sem = std(diff_ts, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(diff_ts),1));
        y_all = [y_all, m+sem, m-sem]; %#ok<AGROW>

        patch([t, fliplr(t)], [m+sem, fliplr(m-sem)], colors.(sgn), 'EdgeColor','none', 'FaceAlpha',0.25);
        plot(t, m, '-', 'Color', colors.(sgn), 'LineWidth', 1.8, 'DisplayName', sprintf('%s-cluster channels', sgn));
    end

    yl = [min(y_all) max(y_all)];
    if isempty(yl) || any(isnan(yl)), yl = [-1 1]; end
    for signtype = {'neg','pos'}
        sgn = signtype{1};
        tr = cluster_time_ranges.(sgn);
        if isempty(tr), continue; end
        patch([tr(1) tr(2) tr(2) tr(1)], [yl(1) yl(1) yl(2) yl(2)], colors.(sgn), ...
            'FaceAlpha', 0.08, 'EdgeColor', colors.(sgn), 'LineStyle', ':', 'HandleVisibility','off');
    end
    yline(0, 'k:', 'HandleVisibility','off');
    xline(3.5, 'k--', 'HandleVisibility','off'); xline(4.0, 'k--', 'HandleVisibility','off');
    xlabel('Time (s, relative to cue onset)'); ylabel('Valid - Invalid alpha power (dB)');
    title(sprintf('%s: neg vs pos cluster time-course (dashed = encoding window 3.5-4.0s)', cue_label), 'FontSize', 10);
    legend('Location','best'); xlim([3.0 4.5]);
    saveas(f, outfile); close(f);
end
