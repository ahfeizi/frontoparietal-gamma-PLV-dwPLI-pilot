% s12_evoked_induced_decomposition.m  --  DIAGNOSTIC. Decomposes the
% encoding-window valid-invalid alpha effect into EVOKED (phase-locked,
% i.e. what a plain ERP would produce when run through a wavelet
% transform) and INDUCED (non-phase-locked oscillatory) components.
%
% WHY: s11's timeline plot showed the neg- and pos-clusters are really
% one continuous, narrow (~250-350ms), sharply biphasic transient
% starting right at Targets onset (3.5s) -- a shape much more typical of
% an evoked visual response (e.g. P1/N1-like) than a sustained
% oscillatory desynchronization. This does NOT mean the effect is an
% artifact (a validity-modulated early evoked response would itself be a
% real and interesting attentional finding, consistent with classic
% Posner-cueing ERP amplitude effects) -- but it changes how the effect
% should be described, and needs to be checked directly rather than
% inferred from shape alone.
%
% METHOD (standard evoked/induced decomposition):
%   TOTAL power(t,f)   = trial-average of |wavelet(single trial)|^2
%   EVOKED power(t,f)  = |wavelet(trial-averaged ERP)|^2
%   INDUCED power(t,f) = TOTAL power(t,f) - EVOKED power(t,f)   [linear units, pre-baseline]
% All three are then dB-baseline-corrected the same way (paper's -500 to
% -100ms pre-cue window) using each measure's own pre-cue values, and the
% valid-invalid DIFFERENCE is computed for each of the three separately.
% If the "evoked" difference curve reproduces the total difference curve
% and "induced" does not, the effect is primarily phase-locked/ERP-driven.
%
% SCOPE: restricted to the alpha-adjacent band (6-14 Hz, for wavelet edge
% margin around the 8-12 Hz band of interest) and to the union of the
% neg+pos cluster channels identified in s09/s11, to keep this tractable
% as a targeted diagnostic rather than a full whole-scalp recomputation.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;

out_dir = fullfile(cfg.deriv_dir, 'alpha_diagnostics');
if ~exist(out_dir, 'dir'); mkdir(out_dir); end

S = load(fullfile(cfg.deriv_dir, 'alpha_group_stats.mat'));

M = load(fullfile(cfg.deriv_dir, 'manifest.mat')); manifest = M.manifest;
srate = cfg.resample_hz;

% ---- restricted analysis params (targeted, not the full s08 sweep) ----
foi_local  = 6:1:14;               % alpha-adjacent band, wavelet edge margin
% toi MUST include the pre-cue baseline window (-0.5 to -0.1s), not just
% the context window around encoding -- otherwise the dB baseline
% correction below finds no baseline samples and silently produces NaN
% everywhere (this bug was caught in an earlier run: the toi previously
% started at 3.0s, entirely excluding the baseline).
toi_local  = -0.5:cfg.alpha.toi_step_s:4.5;  % continuous, covers baseline through context window
bl         = cfg.alpha.baseline_s;
win        = round([-0.5 4.5] * srate);      % epoch samples rel. to cue onset -- covers baseline through context window
lock_col   = cfg.ti.col_cue_sample;

cue_cases = struct('label', {'stick','avatar'});
% ---- cluster channels (union of neg+pos) per cue type, from s09/s11 ----
for ci = 1:2
    stat = eval(sprintf('S.stat_encoding_%s', cue_cases(ci).label));
    chans = {};
    for signtype = {'neg','pos'}
        sgn = signtype{1};
        cf = [sgn 'clusters']; lf = [sgn 'clusterslabelmat'];
        if isfield(stat, cf) && ~isempty(stat.(cf))
            mask = stat.(lf) == 1;
            chan_idx = any(any(mask,2),3);
            chans = union(chans, stat.label(chan_idx));
        end
    end
    cue_cases(ci).chans = chans;
    fprintf('%s: %d union cluster channels for evoked/induced analysis.\n', cue_cases(ci).label, numel(chans));
end

ei_manifest_file = fullfile(cfg.deriv_dir, 'ei_manifest.mat');
done_keys = {};
if exist(ei_manifest_file, 'file')
    dm = load(ei_manifest_file); done_keys = dm.done_keys;
end

for r = 1:height(manifest)
    key  = char(manifest.key(r));
    subj = char(manifest.subject(r));
    out_file = fullfile(cfg.deriv_dir, [key '_evokedInduced.mat']);

    if ismember(key, done_keys) && exist(out_file, 'file')
        fprintf('%s already processed (evoked/induced), skipping (resume).\n', key);
        continue;
    end
    fprintf('\n=== %s (evoked/induced decomposition) ===\n', key);

    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    if isKey(cfg.manual_bad_channels, subj)
        EEG = pop_select(EEG, 'nochannel', cfg.manual_bad_channels(subj));
    end

    ti_path = cfg.authors_deriv_file(subj);
    if ~exist(ti_path, 'file')
        warning('%s: trialinfo not found, skipping.', subj);
        continue;
    end
    TI = load(ti_path, 'trialinfo'); trialinfo = TI.trialinfo;
    nTrials = size(trialinfo, 1);

    % ---- identical inclusion filtering to s08_alpha_epoch_tfr.m ----
    lock = trialinfo(:, lock_col);
    rt_loc  = (trialinfo(:, cfg.ti.col_locresp)  - trialinfo(:, cfg.ti.col_locprobe)) / srate;
    rt_stat = (trialinfo(:, cfg.ti.col_statresp) - trialinfo(:, cfg.ti.col_statq))    / srate;
    mad_loc  = median(abs(rt_loc(~isnan(rt_loc))   - median(rt_loc,'omitnan')));
    mad_stat = median(abs(rt_stat(~isnan(rt_stat)) - median(rt_stat,'omitnan')));
    loc_degenerate  = mad_loc  < cfg.rt.degenerate_mad_thresh_s;
    stat_degenerate = mad_stat < cfg.rt.degenerate_mad_thresh_s;
    rt_ok = true(nTrials,1);
    if cfg.rt.use_rt_filter
        if ~loc_degenerate,  rt_ok = rt_ok & local_rt_ok(rt_loc, cfg);  end
        if ~stat_degenerate, rt_ok = rt_ok & local_rt_ok(rt_stat, cfg); end
    end
    acc_ok = true(nTrials,1);
    if cfg.correct_only
        acc_ok = (trialinfo(:, cfg.accuracy_column) == cfg.accuracy_correct_value);
    end
    boundary_samples = [];
    if isfield(EEG,'event') && ~isempty(EEG.event)
        boundary_samples = [EEG.event(strcmp({EEG.event.type},'boundary')).latency];
    end
    a = lock + win(1); z = lock + win(2);
    in_bounds = (a>=1) & (z<=EEG.pnts);
    crosses_boundary = false(nTrials,1);
    if ~isempty(boundary_samples)
        for t = 1:nTrials
            if in_bounds(t)
                crosses_boundary(t) = any(boundary_samples > a(t) & boundary_samples < z(t));
            end
        end
    end
    good = acc_ok & rt_ok & in_bounds & ~crosses_boundary;
    trl = [a(good), z(good), repmat(win(1), sum(good), 1)];
    ti_good = trialinfo(good, :);

    ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
    cfg_ep = []; cfg_ep.trl = trl;
    epoched = ft_redefinetrial(cfg_ep, ft_raw);
    epoched.trialinfo = ti_good;

    subj_result = struct('subject', subj, 'key', key);

    for ci = 1:numel(cue_cases)
        cue_label = cue_cases(ci).label;
        chans = cue_cases(ci).chans;
        avail_chans = intersect(chans, epoched.label);
        if isempty(avail_chans)
            warning('%s/%s: none of the cluster channels present, skipping.', key, cue_label);
            continue
        end

        cue_val = cfg.ti.(sprintf('cue_%s', cue_label));
        for validity = {'valid','invalid'}
            vlabel = validity{1};
            val_val = cfg.ti.(sprintf('val_%s', vlabel));
            idx = ti_good(:, cfg.ti.col_cue)==cue_val & ti_good(:, cfg.ti.col_val)==val_val;
            if sum(idx) < 5
                warning('%s/%s/%s: <5 trials, skipping this cell.', key, cue_label, vlabel);
                continue
            end

            cfg_sel = []; cfg_sel.trials = find(idx); cfg_sel.channel = avail_chans;
            cond_data = ft_selectdata(cfg_sel, epoched);

            % ---- TOTAL power: trial-average of single-trial power ----
            cfg_tf = []; cfg_tf.method='wavelet'; cfg_tf.output='pow';
            cfg_tf.foi=foi_local; cfg_tf.width=cfg.alpha.cycles; cfg_tf.toi=toi_local;
            cfg_tf.pad='nextpow2'; cfg_tf.keeptrials='no';
            total_tfr = ft_freqanalysis(cfg_tf, cond_data);

            % ---- EVOKED power: power of the trial-averaged ERP ----
            cfg_avg = []; erp = ft_timelockanalysis(cfg_avg, cond_data);
            erp_as_trial = [];
            erp_as_trial.label = erp.label;
            erp_as_trial.time  = {erp.time};
            erp_as_trial.trial = {erp.avg};
            erp_as_trial.fsample = srate;
            erp_as_trial.dimord = 'chan_time';
            evoked_tfr = ft_freqanalysis(cfg_tf, erp_as_trial);

            % ---- align time axes (wavelet edge NaNs can differ slightly) ----
            [~, ia, ib] = intersect(round(total_tfr.time,3), round(evoked_tfr.time,3));
            t_common = total_tfr.time(ia);
            total_pow  = total_tfr.powspctrm(:,:,ia);   % chan x freq x time
            evoked_pow = evoked_tfr.powspctrm(:,:,ib);
            induced_pow = total_pow - evoked_pow;
            induced_pow(induced_pow < 0) = NaN; % negative induced power is a numerical artifact of subtraction, not physical

            bl_idx = t_common >= bl(1) & t_common <= bl(2);
            % CRITICAL FIX: total, evoked, and induced must all be
            % normalized against the SAME baseline (total power's own
            % pre-cue baseline), not each one's own baseline separately.
            % Evoked power's own baseline is much smaller than total's
            % (ERP averaging suppresses non-phase-locked noise, so the
            % pre-cue "evoked" baseline is close to the noise floor),
            % so self-normalizing evoked into dB inflates it artificially
            % relative to total -- this was caught by an impossible
            % evoked/total ratio >100% (evoked can never exceed total,
            % by Jensen's inequality: |mean(x)|^2 <= mean(|x|^2)).
            %
            % Using a single common baseline (total's) means the three
            % quantities, expressed as %-of-baseline-power, satisfy
            % total_pct = evoked_pct + induced_pct EXACTLY (linear
            % normalization preserves the total=evoked+induced identity;
            % dB/log does not). dB is kept only for the 'total' line
            % (matching the rest of the pipeline's convention); evoked
            % and induced are reported as % change from total's baseline.
            total_bl = mean(total_pow(:,:,bl_idx), 3, 'omitnan');  % chan x freq

            total_db    = 10*log10(total_pow ./ total_bl);                  % dB, standard
            evoked_pct  = 100 * (evoked_pow  ./ total_bl - 1);              % % change from total's baseline
            induced_pct = 100 * (induced_pow ./ total_bl - 1);
            total_pct   = 100 * (total_pow   ./ total_bl - 1);              % same info as total_db, linear units

            % trim back down to the context window for storage -- the
            % pre-cue span was only needed transiently for baseline calc
            ctx_idx = t_common >= 3.0;
            t_ctx = t_common(ctx_idx);
            total_db    = total_db(:,:,ctx_idx);
            total_pct   = total_pct(:,:,ctx_idx);
            evoked_pct  = evoked_pct(:,:,ctx_idx);
            induced_pct = induced_pct(:,:,ctx_idx);

            fname = sprintf('%s_%s', cue_label, vlabel);
            subj_result.(fname) = struct('time', t_ctx, 'freq', total_tfr.freq, 'chans', {avail_chans}, ...
                'total_db', total_db, 'total_pct', total_pct, 'evoked_pct', evoked_pct, 'induced_pct', induced_pct);
        end
    end

    save(out_file, 'subj_result');
    done_keys{end+1} = key; %#ok<AGROW>
    save(ei_manifest_file, 'done_keys');
end

fprintf('\nPer-subject evoked/induced decomposition complete. Aggregating...\n');

%% ---- group-level aggregation + plot ----
% Collected as a CELL array, not a struct array: subjects can have
% different subsets of fields present (e.g. a cell skipped for <5 trials
% or missing cluster channels), and concatenating structs with unequal
% fields errors ("dissimilar structures"). A cell array has no such
% requirement.
all_ei = {};
for r = 1:height(manifest)
    key = char(manifest.key(r));
    f = fullfile(cfg.deriv_dir, [key '_evokedInduced.mat']);
    if ~exist(f, 'file'), continue; end
    L = load(f);
    all_ei{end+1} = L.subj_result; %#ok<AGROW>
end
fprintf('Loaded evoked/induced results for %d subjects.\n', numel(all_ei));

report_lines = {};
report_lines{end+1} = sprintf('Evoked/induced decomposition: %s', datestr(now));

for ci = 1:numel(cue_cases)
    cue_label = cue_cases(ci).label;
    vfield = sprintf('%s_valid', cue_label);
    ifield = sprintf('%s_invalid', cue_label);

    has_both = cellfun(@(x) isfield(x, vfield) && isfield(x, ifield) && ~isempty(x.(vfield)) && ~isempty(x.(ifield)), all_ei);
    subs = all_ei(has_both);
    if isempty(subs)
        fprintf('%s: no subjects with both conditions present, skipping.\n', cue_label);
        continue
    end

    t = subs{1}.(vfield).time;
    freq_mask = subs{1}.(vfield).freq >= 8 & subs{1}.(vfield).freq <= 12;

    diff_total = nan(numel(subs), numel(t));
    diff_evoked = nan(numel(subs), numel(t));
    diff_induced = nan(numel(subs), numel(t));
    for s = 1:numel(subs)
        v = subs{s}.(vfield); iv = subs{s}.(ifield);
        % all three use the % change (from total's own baseline) fields,
        % so total_diff = evoked_diff + induced_diff holds exactly
        % (unlike dB, which is not additive) -- this is the actual
        % evoked/induced decomposition of the effect, in comparable units.
        diff_total(s,:)   = squeeze(mean(v.total_pct(:,freq_mask,:)   - iv.total_pct(:,freq_mask,:),   [1 2], 'omitnan'));
        diff_evoked(s,:)  = squeeze(mean(v.evoked_pct(:,freq_mask,:)  - iv.evoked_pct(:,freq_mask,:),  [1 2], 'omitnan'));
        diff_induced(s,:) = squeeze(mean(v.induced_pct(:,freq_mask,:) - iv.induced_pct(:,freq_mask,:), [1 2], 'omitnan'));
    end

    plot_three_way(t, diff_total, diff_evoked, diff_induced, cue_label, ...
        fullfile(out_dir, sprintf('evoked_induced_%s.png', cue_label)));

    win_idx = t >= 3.5 & t <= 4.0;
    peak_total   = max(abs(mean(diff_total(:,win_idx), 1, 'omitnan')));
    peak_evoked  = max(abs(mean(diff_evoked(:,win_idx), 1, 'omitnan')));
    peak_induced = max(abs(mean(diff_induced(:,win_idx), 1, 'omitnan')));
    fprintf('%s: peak |valid-invalid| in encoding window (%% of total baseline) -- total=%.2f%%, evoked=%.2f%%, induced=%.2f%% (evoked share=%.0f%%)\n', ...
        cue_label, peak_total, peak_evoked, peak_induced, 100*peak_evoked/peak_total);
    report_lines{end+1} = sprintf('%s: peak diff (%% of total baseline) -- total=%.2f evoked=%.2f induced=%.2f (evoked share=%.0f%%)', ...
        cue_label, peak_total, peak_evoked, peak_induced, 100*peak_evoked/peak_total); %#ok<AGROW>
end

report_file = fullfile(out_dir, 'evoked_induced_report.txt');
fid = fopen(report_file,'w'); fprintf(fid, '%s\n', report_lines{:}); fclose(fid);
fprintf('\nSaved figures + report to %s\n', out_dir);

% ======================= local functions =======================
function ok = local_rt_ok(rt, cfg)
    ok = ~isnan(rt) & (rt >= cfg.rt.min_rt_s);
    valid_rt = rt(ok);
    if numel(valid_rt) >= 5
        med = median(valid_rt); madv = median(abs(valid_rt - med));
        if madv > 0
            lo = med - cfg.rt.outlier_n*madv*1.4826; hi = med + cfg.rt.outlier_n*madv*1.4826;
            ok = ok & (rt >= lo) & (rt <= hi);
        end
    end
end

function plot_three_way(t, diff_total, diff_evoked, diff_induced, cue_label, outfile)
    f = figure('Visible','off','Color','w','Position',[100 100 750 420]); hold on;
    specs = {diff_total,'total',[0 0 0]; diff_evoked,'evoked',[0.85 0.2 0.1]; diff_induced,'induced',[0.1 0.4 0.85]};
    for i = 1:3
        d = specs{i,1}; m = mean(d,1,'omitnan'); sem = std(d,0,1,'omitnan')./sqrt(sum(~isnan(d),1));
        patch([t fliplr(t)], [m+sem fliplr(m-sem)], specs{i,3}, 'EdgeColor','none','FaceAlpha',0.15,'HandleVisibility','off');
        plot(t, m, '-', 'Color', specs{i,3}, 'LineWidth', 1.8, 'DisplayName', specs{i,2});
    end
    yline(0,'k:','HandleVisibility','off');
    xline(3.5,'k--','HandleVisibility','off'); xline(4.0,'k--','HandleVisibility','off');
    xlabel('Time (s, relative to cue onset)'); ylabel('Valid - Invalid (% of total pre-cue baseline power)');
    title(sprintf('%s: total = evoked + induced (dashed = encoding window)', cue_label), 'FontSize', 10);
    legend('Location','best'); xlim([t(1) t(end)]);
    saveas(f, outfile); close(f);
end
