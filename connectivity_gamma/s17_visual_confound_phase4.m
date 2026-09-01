% s17_visual_confound_phase4.m -- Phase 4: does the avatar vs stick cue
% differ in low-level visual response (occipital ERP amplitude) during the
% shift transition, independent of any cognitive/attentional processing?
%
% RATIONALE: avatar and stick cues differ physically (shape, articulated
% motion vs rigid rotation), so any d_cue effect (main effect of cue type)
% could in principle reflect this low-level visual difference rather than
% social/non-social processing per se. This does NOT threaten the
% preprint's central finding (validity effect, d_val) unless it interacts
% with validity -- and d_int has consistently been null/borderline-null
% across all checks so far (S10, RT_stat sensitivity). This check is about
% d_cue and general robustness/transparency, not about undermining d_val.
%
% SHIFT-WINDOW DERIVATION (inferred, not a direct trialinfo column):
% trialinfo has Cue onset (col 9, t=0 reference) and Targets onset (col 10,
% ~3.5s post-cue per the task paper's fixed timeline), but no separate
% "shift onset" column. Per Gregory, Wang & Kessler (2022) Fig. 1 / Methods,
% the cue shifts left/right during the 500 ms immediately before targets
% appear (shift: 3.0-3.5s post cue onset, targets at 3.5s). So:
%   shift_window = [Targets_sample - round(0.5*srate), Targets_sample - 1]
% SANITY CHECK (printed before trusting this): median/range of
% (Targets_sample - Cue_sample)/srate across all trials, which should be
% close to 3.5s if the task's fixed timeline holds for this dataset.
%
% ROI: occipital channels {O1, O2, Oz, POz} -- closest to primary visual
% cortex, standard choice for a low-level visual (not cognitive) proxy.
%
% Trials: ALL trials (not just "good"), since a low-level visual response
% should not depend on later accuracy/RT -- restricting would bias toward
% a specific outcome subgroup for no principled reason.
%
% Analysis: paired comparison (avatar vs stick) of occipital broadband ERP
% amplitude (RMS in the shift window, baseline-corrected) per subject, plus
% a cluster-based permutation time-course as a secondary/visual check.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

occ_channels = {'O1','O2','Oz','POz'};
srate = cfg.resample_hz;
baseline_s = [-0.2 0];   % relative to shift onset

rms_rows = {};
timecourse_avatar = {}; timecourse_stick = {};   % for optional cluster permutation
gap_check = [];

for su = 1:height(manifest)
    subj = char(manifest.subject(su)); key = char(manifest.key(su));
    ti_path = cfg.authors_deriv_file(subj);
    if ~isfile(ti_path)
        warning('%s: trialinfo not found, skipping.', subj); continue;
    end
    TI = load(ti_path, 'trialinfo'); trialinfo = TI.trialinfo;

    cue_sample = trialinfo(:, cfg.ti.col_cue_sample);
    targets_sample = trialinfo(:, cfg.ti.col_targets_sample);
    cueV = trialinfo(:, cfg.ti.col_cue);

    gap_s = (targets_sample - cue_sample) / srate;
    gap_check = [gap_check; gap_s]; %#ok<AGROW>

    set_file = fullfile(cfg.deriv_dir, [key '_clean.set']);
    if ~isfile(set_file), warning('%s: clean.set not found, skipping.', subj); continue; end
    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    occIdx = util_chanidx(EEG, occ_channels);
    if isempty(occIdx), warning('%s: no occipital ROI channels found, skipping.', subj); continue; end

    shift_len = round(0.5*srate);
    base_len_pre = round(abs(baseline_s(1))*srate);

    for c = 1:2
        cueLab = {'stick','avatar'}; cv = [cfg.ti.cue_stick, cfg.ti.cue_avatar];
        sel = find(cueV == cv(c));
        amp_trials = nan(numel(sel),1);
        erp_accum = zeros(1, shift_len);
        n_erp = 0;
        for t = 1:numel(sel)
            ts = targets_sample(sel(t));
            shift_start = ts - shift_len;
            base_start = shift_start - base_len_pre;
            if base_start < 1 || ts > EEG.pnts, continue; end

            baseline_seg = mean(EEG.data(occIdx, base_start:shift_start-1), 1);
            baseline_mean = mean(baseline_seg);
            shift_seg = mean(EEG.data(occIdx, shift_start:ts-1), 1) - baseline_mean;

            amp_trials(t) = rms(shift_seg);
            erp_accum = erp_accum + shift_seg; n_erp = n_erp + 1;
        end
        rms_rows(end+1,:) = {subj, cueLab{c}, mean(amp_trials,'omitnan')}; %#ok<SAGROW>
        if n_erp > 0
            if c==1, timecourse_stick{end+1} = erp_accum/n_erp; %#ok<SAGROW>
            else,    timecourse_avatar{end+1} = erp_accum/n_erp; %#ok<SAGROW>
            end
        end
    end
    fprintf('  %s done\n', subj);
end

fprintf('\nSANITY CHECK -- Targets-Cue gap: median=%.3fs, range=[%.3f %.3f]s (expect ~3.5s if task timeline is fixed)\n', ...
    median(gap_check,'omitnan'), min(gap_check), max(gap_check));
if abs(median(gap_check,'omitnan') - 3.5) > 0.1
    warning(['Targets-Cue gap median deviates from the assumed 3.5s by more than 100ms -- ' ...
        'the shift_window derivation may not hold for this dataset. Verify before trusting results below.']);
end

% ======================= paired comparison: occipital RMS amplitude =======================
Trms = cell2table(rms_rows, 'VariableNames', {'subject','cue','rms_amp'});
subs = unique(Trms.subject, 'stable');
av = nan(numel(subs),1); st = nan(numel(subs),1);
for i=1:numel(subs)
    ra = Trms(strcmp(Trms.subject,subs{i}) & strcmp(Trms.cue,'avatar'),:);
    rs = Trms(strcmp(Trms.subject,subs{i}) & strcmp(Trms.cue,'stick'),:);
    if ~isempty(ra), av(i)=ra.rms_amp(1); end
    if ~isempty(rs), st(i)=rs.rms_amp(1); end
end
d = av - st; d = d(~isnan(d)); n = numel(d);
dz = mean(d)/std(d); t_obs = mean(d)/(std(d)/sqrt(n)); p = 2*(1-tcdf(abs(t_obs), n-1));

fprintf('\n--- Occipital broadband RMS amplitude, avatar vs stick (shift window) ---\n');
fprintf('n=%d, avatar mean=%.3f, stick mean=%.3f, dz=%.3f, t=%.3f, p=%.4f\n', ...
    n, mean(av,'omitnan'), mean(st,'omitnan'), dz, t_obs, p);
if p < 0.05
    fprintf(['INTERPRETATION: significant low-level visual amplitude difference between cue types.\n' ...
        'This is expected (avatar and stick are physically different stimuli) and does not threaten\n' ...
        'the validity (d_val) finding unless it interacts with validity -- d_int has been consistently\n' ...
        'null/borderline across all checks. Report as an acknowledged, expected visual difference,\n' ...
        'not a confound of the reported effect.\n']);
else
    fprintf('INTERPRETATION: no significant difference detected in this broadband amplitude proxy.\n');
end

writetable(Trms, fullfile(cfg.deriv_dir, 'visual_confound_occipital_rms.csv'));
fprintf('Saved visual_confound_occipital_rms.csv\n');

% ======================= secondary: cluster-based permutation time-course =======================
if ~isempty(timecourse_avatar) && ~isempty(timecourse_stick)
    n2 = min(numel(timecourse_avatar), numel(timecourse_stick));
    Ma = cell2mat(timecourse_avatar(1:n2)'); Ms = cell2mat(timecourse_stick(1:n2)');
    Dtc = Ma - Ms;   % [subject x time]
    tobs_tc = mean(Dtc) ./ (std(Dtc)./sqrt(n2));
    rng(1); nperm = cfg.stat.n_perm;
    alpha_form = 0.05; t_thr = tinv(1-alpha_form/2, n2-1);
    maxmass_null = zeros(nperm,1);
    for p_i = 1:nperm
        sgn = sign(rand(n2,1)-0.5); Dp = Dtc.*sgn; tp = mean(Dp)./(std(Dp)./sqrt(n2));
        supra = abs(tp) > t_thr;
        if any(supra)
            % simple contiguous-run cluster mass (1D time series)
            d_supra = diff([0 supra 0]); starts = find(d_supra==1); ends = find(d_supra==-1)-1;
            masses = arrayfun(@(a,b) sum(abs(tp(a:b))), starts, ends);
            maxmass_null(p_i) = max(masses);
        end
    end
    supra_obs = abs(tobs_tc) > t_thr;
    d_supra = diff([0 supra_obs 0]); starts = find(d_supra==1); ends = find(d_supra==-1)-1;
    fprintf('\n--- Time-course cluster permutation (occipital ERP, shift window, n=%d) ---\n', n2);
    if isempty(starts)
        fprintf('No supra-threshold time cluster found.\n');
    else
        for k = 1:numel(starts)
            m = sum(abs(tobs_tc(starts(k):ends(k))));
            p_clust = mean(maxmass_null >= m);
            fprintf('  cluster %d: samples %d-%d (%.0f-%.0fms), mass=%.2f, p_cluster=%.4f\n', ...
                k, starts(k), ends(k), (starts(k)-1)/srate*1000, (ends(k)-1)/srate*1000, m, p_clust);
        end
    end
    save(fullfile(cfg.deriv_dir, 'visual_confound_timecourse.mat'), 'Ma','Ms','tobs_tc','srate');
    fprintf('Saved visual_confound_timecourse.mat\n');
else
    fprintf('\nInsufficient data for time-course cluster permutation.\n');
end
