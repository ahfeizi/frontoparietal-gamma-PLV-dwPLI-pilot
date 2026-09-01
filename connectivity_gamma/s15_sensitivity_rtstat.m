% s15_sensitivity_rtstat.m -- Sensitivity check: the RT_stat filter
% (cfg.rt.use_rt_filter) was auto-disabled for 43/47 subjects because their
% RT_stat was degenerate (near-constant, logged in rt_qc_flags.csv). It was
% only ACTUALLY applied for the remaining 4 subjects -- confirmed from
% rt_qc_flags.csv: sub-01, sub-02, sub-03, sub-04 (rt_stat_degenerate=0);
% all others =1. This creates an inconsistent trial-inclusion criterion
% across subjects. This script tests whether disabling the RT_stat filter
% for those same 4 subjects (making trial inclusion uniform across all 47)
% changes the key gamma dwPLI/PLV results.
%
% RECONSTRUCTION OF "GOOD" TRIALS UNDER THE ALTERNATE POLICY -- derived
% analytically from trials.mat alone (no need for raw trialinfo/boundary
% info again), using the fact that s02_epoch.m's 'reason' field records
% the FIRST failing criterion in a fixed priority order (out_of_bounds >
% crosses_boundary > incorrect > rt_outlier_or_anticipatory):
%   - reason in {'out_of_bounds','crosses_boundary'}: exclusion is
%     unrelated to RT -> good2 = false, unchanged.
%   - reason == 'incorrect': exclusion is accuracy-based, unrelated to RT
%     -> good2 = false, unchanged.
%   - reason == 'rt_outlier_or_anticipatory' OR good==true: this trial
%     passed in_bounds & ~crosses_boundary & acc_ok, so RT is the only
%     thing determining good2. Recompute rt_ok_loc_only using the EXACT
%     same local_rt_ok() logic as s02 (copied verbatim), applied ONLY to
%     rt_loc (dropping the rt_stat term) -> good2 = rt_ok_loc_only(t).
% This is an exact reconstruction, not an approximation, because
% local_rt_ok() is a pure function of one subject's rt array.
%
% Connectivity (dwPLI/PLV, gamma_high only -- where the key effect lives)
% is then RECOMPUTED for these 4 subjects only, using the same recipe as
% s03_connectivity.m (bandpass + Hilbert, NOT wavelet), and merged with the
% other 43 subjects' unchanged values from conn_results.mat.

clear; clc;
cfg = config();
rng(1);
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;
RTflags = readtable(fullfile(cfg.deriv_dir,'rt_qc_flags.csv'));

affected_subs = RTflags.subject(RTflags.rt_stat_degenerate == 0);
fprintf('Subjects where RT_stat filter was actually applied (to be re-run without it): %s\n', ...
    strjoin(affected_subs, ', '));

bn = 'gamma_high';
band_range = cfg.bands.(bn);
nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_parietal);
win = round(cfg.maint_window * cfg.resample_hz);
wlen = win(2) - win(1) + 1;
cueLabs = {'stick','avatar'}; cueVals = [cfg.ti.cue_stick, cfg.ti.cue_avatar];
valLabs = {'valid','invalid'}; valVals = [cfg.ti.val_valid, cfg.ti.val_invalid];

R2 = R;   % copy; we will overwrite rows for affected_subs, band==gamma_high

for as = 1:numel(affected_subs)
    subj = affected_subs{as};
    key = subj;   % ASSUMPTION confirmed working in s13 -- file key == subject id
    trials_file = fullfile(cfg.deriv_dir, [key '_trials.mat']);
    T = load(trials_file); trials = T.trials;

    rt_loc_all = [trials.rt_loc];
    rt_ok_loc_only = local_rt_ok(rt_loc_all, cfg);

    good2 = false(1, numel(trials));
    for t = 1:numel(trials)
        rsn = trials(t).reason;
        if trials(t).good || strcmp(rsn, 'rt_outlier_or_anticipatory')
            good2(t) = rt_ok_loc_only(t);
        else
            good2(t) = false;   % out_of_bounds / crosses_boundary / incorrect -- unchanged
        end
    end
    n_added = sum(good2) - sum([trials.good]);
    fprintf('%s: good (orig)=%d, good2 (RT_stat-disabled)=%d (delta=%+d)\n', ...
        subj, sum([trials.good]), sum(good2), n_added);

    % ---- recompute dwPLI/PLV for gamma_high under good2, this subject only ----
    set_file = fullfile(cfg.deriv_dir, sprintf('%s_clean.set', key));
    EEG = pop_loadset('filename', set_file);
    fr = util_chanidx(EEG, cfg.roi_frontal);
    pa = util_chanidx(EEG, cfg.roi_parietal);

    Xf = util_bandpass(double(EEG.data), EEG.srate, band_range, cfg.filt_order_factor);
    Z  = hilbert(Xf.').';

    good2_trials = trials(good2);
    for c = 1:numel(cueLabs)
        for v = 1:numel(valLabs)
            sel = good2_trials([good2_trials.cue]==cueVals(c) & [good2_trials.valid]==valVals(v));
            ntr = numel(sel);
            if ntr < 5
                fprintf('  %s-%s: only %d good2 trials, skipping cell (leaving conn_results.mat value in place)\n', ...
                    cueLabs{c}, valLabs{v}, ntr);
                continue;
            end
            seg = nan(EEG.nbchan, wlen, ntr); kept = 0;
            for t = 1:ntr
                a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                if a < 1 || z > size(Z,2), continue; end
                kept = kept + 1; seg(:,:,kept) = Z(:, a:z);
            end
            seg = seg(:,:,1:kept);
            if kept < 5, continue; end

            np = numel(fr)*numel(pa);
            dwpli_pairs = zeros(1,np); plv_pairs = zeros(1,np); p = 0;
            for i = 1:numel(fr)
                zi = reshape(seg(fr(i),:,:), wlen, kept);
                for j = 1:numel(pa)
                    p = p+1;
                    zj = reshape(seg(pa(j),:,:), wlen, kept);
                    X = zi.*conj(zj); im = imag(X);
                    s_im=sum(im,2); s_abs=sum(abs(im),2); s_sq=sum(im.^2,2);
                    dwpli_t = (s_im.^2 - s_sq) ./ (s_abs.^2 - s_sq);
                    plv_t = abs(mean(X./abs(X), 2));
                    dwpli_pairs(p) = mean(dwpli_t,'omitnan');
                    plv_pairs(p) = mean(plv_t,'omitnan');
                end
            end

            row_mask = strcmp(R2.subject,subj) & strcmp(R2.band,bn) & strcmp(R2.cue,cueLabs{c}) & strcmp(R2.valid,valLabs{v});
            if sum(row_mask) ~= 1
                warning('%s %s-%s: expected exactly 1 matching row in conn_results, found %d -- skipping update.', subj, cueLabs{c}, valLabs{v}, sum(row_mask));
                continue;
            end
            R2.dwpli(row_mask) = mean(dwpli_pairs,'omitnan');
            R2.plv(row_mask) = mean(plv_pairs,'omitnan');
            R2.dwpli_pairs(row_mask) = {dwpli_pairs};
            R2.plv_pairs(row_mask) = {plv_pairs};
            R2.ntrials(row_mask) = kept;
        end
    end
    fprintf('  %s recomputed for gamma_high.\n', subj);
end

% ======================= compare original vs RT_stat-disabled table (gamma_high only) =======================
Delta = 0.4; r_prior = 0.707;
metrics = {'dwpli','plv'};
getval = @(Ruse,subj,cueLab,valLab,mt) local_getval(Ruse, subj, bn, mt, cueLab, valLab);

rows_out = {};   % <-- NEW: accumulator for the comparison table

for which = 1:2
    if which==1, Ruse = R;  tag = 'ORIGINAL (RT_stat filter as pipeline ran it)';
    else,        Ruse = R2; tag = 'RT_stat filter DISABLED for sub-01..04';
    end
    subs = unique(Ruse.subject);
    fprintf('\n--- %s ---\n', tag);
    fprintf('%-8s %-8s %3s  %7s  %7s  %8s  %8s  %8s\n','contrast','metric','n','dz','t','BF10','BF01','p_TOST');
    for m = 1:numel(metrics)
        mt = metrics{m};
        av=nan(numel(subs),1); ai=nan(numel(subs),1); sv=nan(numel(subs),1); si_=nan(numel(subs),1);
        for i=1:numel(subs)
            av(i)=getval(Ruse,subs{i},'avatar','valid',mt);
            ai(i)=getval(Ruse,subs{i},'avatar','invalid',mt);
            sv(i)=getval(Ruse,subs{i},'stick','valid',mt);
            si_(i)=getval(Ruse,subs{i},'stick','invalid',mt);
        end
        contrasts = struct('d_cue',(av+ai)/2-(sv+si_)/2, 'd_val',(av+sv)/2-(ai+si_)/2, 'd_int',(av-ai)-(sv-si_));
        cn = fieldnames(contrasts);
        for c = 1:numel(cn)
            d = contrasts.(cn{c}); d = d(~isnan(d));
            r0 = local_eval(d, Delta, r_prior);
            fprintf('%-8s %-8s %3d  %7.3f  %7.2f  %8.3f  %8.3f  %8.4f\n', cn{c}, mt, r0{1}, r0{2}, r0{3}, r0{5}, r0{6}, r0{7});
            rows_out(end+1,:) = {tag, cn{c}, mt, r0{1}, r0{2}, r0{3}, r0{4}, r0{5}, r0{6}, r0{7}, r0{8}}; %#ok<SAGROW>  % <-- NEW
        end
    end
end

% <-- NEW: save the comparison table
Tout = cell2table(rows_out, 'VariableNames', ...
    {'policy','contrast','metric','n','dz','t','p_freq','BF10','BF01','p_TOST','equivalent_at_0p4'});
writetable(Tout, fullfile(cfg.deriv_dir, 'sensitivity_rtstat_gammahigh.csv'));
fprintf('\nSaved sensitivity_rtstat_gammahigh.csv\n');

save(fullfile(cfg.deriv_dir, 'conn_results_rtstat_sensitivity.mat'), 'R2');
fprintf('\nSaved conn_results_rtstat_sensitivity.mat (R2, with sub-01..04 gamma_high recomputed under RT_stat-disabled policy)\n');

% ======================= local functions =======================
function v = local_getval(R, subj, bn, mt, cueLab, valLab)
    r = R(strcmp(R.subject,subj) & strcmp(R.band,bn) & strcmp(R.cue,cueLab) & strcmp(R.valid,valLab), :);
    if isempty(r), v = NaN; else, v = r.(mt)(1); end
end

function row = local_eval(d, Delta, r_prior)
    n = numel(d); dz = mean(d)/std(d); se = std(d)/sqrt(n); tobs = mean(d)/se; v = n-1;
    p_freq = 2*(1-tcdf(abs(tobs),v));
    BF10 = local_jzs_bf10(tobs, n, r_prior); BF01 = 1/BF10;
    bound_raw = Delta*std(d);
    t_lo=(mean(d)-(-bound_raw))/se; t_hi=(mean(d)-bound_raw)/se;
    p_tost = max(1-tcdf(t_lo,v), tcdf(t_hi,v));
    row = {n, dz, tobs, p_freq, BF10, BF01, p_tost, p_tost<0.05};
end

function bf10 = local_jzs_bf10(t, n, r)
    v = n-1; num = (1+t^2/v)^(-(v+1)/2);
    integrand = @(g) (1+n.*g).^(-0.5).*(1+t^2./((1+n.*g).*v)).^(-(v+1)/2).*(r/sqrt(2*pi)).*g.^(-1.5).*exp(-(r^2)./(2.*g));
    denom = integral(integrand,0,Inf); bf10 = 1/(num/denom);
end

function ok = local_rt_ok(rt, cfg)
    % Copied verbatim from s02_epoch.m -- same MAD-based outlier logic,
    % applied here only to rt_loc (the rt_stat term is intentionally
    % dropped for this sensitivity check).
    ok = ~isnan(rt) & (rt >= cfg.rt.min_rt_s);
    valid_rt = rt(ok);
    if numel(valid_rt) >= 5
        med = median(valid_rt); madv = median(abs(valid_rt-med));
        if madv > 0
            lo = med - cfg.rt.outlier_n*madv*1.4826;
            hi = med + cfg.rt.outlier_n*madv*1.4826;
            ok = ok & (rt>=lo) & (rt<=hi);
        end
    end
end
