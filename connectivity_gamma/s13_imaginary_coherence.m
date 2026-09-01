% s13_imaginary_coherence.m -- Reviewer-motivated check (open item from
% Pilot 2 revision plan): is the valid_main / PLV / gamma_high frontal-
% parietal effect (F2-P4 and the surrounding cluster, s06/s06c/s10) real
% connectivity, or a volume-conduction / common-reference artifact?
%
% RATIONALE: PLV is sensitive to zero-phase-lag (volume conduction, shared
% reference) contributions; dwPLI is constructed specifically to remove
% them, which is exactly why the same contrast was null in dwPLI
% (BF01=6.13, s10) while strong in PLV (BF10=45.16). That dissociation
% ALONE is already suggestive of a volume-conduction contribution, but
% this script provides a direct, independent check: imaginary coherency
% (ImC = imag(Sxy / sqrt(Sxx*Syy))) is, like dwPLI, blind to zero-lag
% coupling. If the validity effect survives in ImC, that argues for a
% genuine (non-zero-lag) synchronization effect specific to PLV's
% sensitivity/normalization rather than to volume conduction per se. If it
% vanishes in ImC too, that strengthens the volume-conduction concern.
%
% RECONSTRUCTION RECIPE -- copied exactly as specified (not re-derived):
% no complex Fourier/Hilbert coefficients were saved by s03_connectivity.m
% (Z and seg were computed in-memory and discarded); this script rebuilds
% them from the same two upstream files s03 uses:
%   <key>_clean.set   (s01_preprocess.m)  -- continuous preprocessed EEG
%   <key>_trials.mat  (s02_epoch.m)       -- trial table (lock_sample, cue,
%                                            valid, good, ...)
% via bandpass filter + Hilbert transform (NOT wavelet -- confirmed no
% Morlet decomposition exists anywhere in the connectivity pipeline; that
% was Pilot 1's s07_grandavg_ersp.m only, unrelated here).
%
% *** ASSUMPTION NEEDING CONFIRMATION -- flagged, not guessed ***
% Subject file-key format is assumed to match R.subject values from
% conn_results.mat (e.g. 'sub-01') used directly as the '<key>' in
% '<key>_clean.set' / '<key>_trials.mat'. If the on-disk filenames differ
% (e.g. zero-padding, a task suffix), fix local_subject_files() below.

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');   % just to get the subject list
subs = unique(R.subject, 'stable');

bn = 'gamma_high';          % the band with the signal; rerun with 'gamma_low' if wanted
band_range = cfg.bands.(bn);

nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_parietal);
labs = cell(1,nF*nT);
for i=1:nF, for j=1:nT, labs{(i-1)*nT+j} = sprintf('%s-%s', cfg.roi_frontal{i}, cfg.roi_parietal{j}); end, end

conditions = {'avatar','valid'; 'avatar','invalid'; 'stick','valid'; 'stick','invalid'};
nCond = size(conditions,1);
% Numeric codes copied verbatim from config.m / s03_connectivity.m --
% trials.cue and trials.valid are numeric, NOT strings:
%   cfg.ti.cue_stick=1, cfg.ti.cue_avatar=2, cfg.ti.val_valid=1, cfg.ti.val_invalid=2
cue_code   = struct('avatar', cfg.ti.cue_avatar, 'stick', cfg.ti.cue_stick);
valid_code = struct('valid',  cfg.ti.val_valid,  'invalid', cfg.ti.val_invalid);

% ImC(subject, condition, pair) and, for the volume-conduction diagnostic,
% the pooled circular mean/resultant length of the trial-level phase
% difference at F2-P4 specifically.
ImC = nan(numel(subs), nCond, nF*nT);
target_pair = find(strcmp(labs,'F2-P4'));
phase_diag = struct('subject',{},'condition',{},'mean_angle_deg',{},'resultant_len',{},'n_trials',{});

fprintf('Reconstructing analytic signal and computing imaginary coherency (%s)...\n', bn);
n_ok = 0; n_skipped_files = 0;
for si = 1:numel(subs)
    key = subs{si};  % ASSUMPTION: file key == R.subject value, see header note
    [set_file, trials_file] = local_subject_files(cfg, key);
    if ~isfile(set_file) || ~isfile(trials_file)
        warning('Missing files for %s (%s / %s) -- skipping subject.', key, set_file, trials_file);
        n_skipped_files = n_skipped_files + 1;
        continue;
    end

    EEG = pop_loadset('filename', set_file);
    Tt = load(trials_file); trials = Tt.trials;

    Xf = util_bandpass(double(EEG.data), EEG.srate, band_range, cfg.filt_order_factor);
    Z  = hilbert(Xf.').';   % [nchan x ntime] continuous analytic signal

    win = round(cfg.maint_window * cfg.resample_hz);
    wlen = win(2) - win(1) + 1;
    good_trials = trials([trials.good]);

    seg = nan(EEG.nbchan, wlen, numel(good_trials));
    n_out_of_range = 0;
    for t = 1:numel(good_trials)
        a = good_trials(t).lock_sample + win(1);
        z = a + wlen - 1;
        if a < 1 || z > size(Z,2), n_out_of_range = n_out_of_range + 1; continue; end   % guard against edge trials
        seg(:,:,t) = Z(:, a:z);
    end
    if n_out_of_range > 0
        fprintf('    %s: %d/%d trials out of range (window exceeds recording bounds)\n', ...
            key, n_out_of_range, numel(good_trials));
    end

    cueArr   = [good_trials.cue];
    validArr = [good_trials.valid];
    if si == 1
        fprintf('    DIAGNOSTIC (sub 1 only) -- unique cue values: %s | unique valid values: %s\n', ...
            mat2str(unique(cueArr)), mat2str(unique(validArr)));
    end

    fIdx = util_chanidx(EEG, cfg.roi_frontal);
    pIdx = util_chanidx(EEG, cfg.roi_parietal);

    for ci = 1:nCond
        cueLab = conditions{ci,1}; valLab = conditions{ci,2};
        sel = (cueArr == cue_code.(cueLab)) & (validArr == valid_code.(valLab));
        segc = seg(:,:,sel);   % [nchan x wlen x ntrial_cond]
        % drop any trial that has NaN anywhere in it (out-of-range window),
        % rather than letting mean() silently propagate NaN through everything
        bad_trial = squeeze(any(any(isnan(segc),1),2));
        segc = segc(:,:,~bad_trial);
        ntr = size(segc,3);
        if si==1
            fprintf('    condition %s/%s: %d trials selected, %d usable after NaN-trial drop\n', ...
                cueLab, valLab, sum(sel), ntr);
        end
        if ntr < 3, continue; end

        for i = 1:nF
            zf = squeeze(segc(fIdx(i),:,:));   % [wlen x ntrial]
            for j = 1:nT
                zp = squeeze(segc(pIdx(j),:,:));
                Sxy = mean(mean(zf .* conj(zp)));
                Sxx = mean(mean(abs(zf).^2));
                Syy = mean(mean(abs(zp).^2));
                coh = Sxy / sqrt(Sxx*Syy);
                k = (i-1)*nT + j;
                ImC(si,ci,k) = imag(coh);

                if k == target_pair
                    dphi = angle(zf) - angle(zp);          % [wlen x ntrial]
                    dphi = dphi(:);
                    mvec = mean(exp(1i*dphi));
                    phase_diag(end+1) = struct('subject',key,'condition',sprintf('%s_%s',cueLab,valLab), ...
                        'mean_angle_deg', angle(mvec)*180/pi, 'resultant_len', abs(mvec), ...
                        'n_trials', ntr); %#ok<SAGROW>
                end
            end
        end
    end
    fprintf('  %s done\n', key);
    n_ok = n_ok + 1;
end
fprintf('\nSubject processing summary: %d ok, %d skipped (missing files), %d total.\n', ...
    n_ok, n_skipped_files, numel(subs));

% ======================= d_val contrast on ImC, ROI-average BF/TOST =======================
Delta = 0.4; r_prior = 0.707;   % copied from s12/s10, unchanged
av = squeeze(ImC(:,1,:)); ai = squeeze(ImC(:,2,:));
sv = squeeze(ImC(:,3,:)); si_ = squeeze(ImC(:,4,:));

roi_d_val = nanmean( (av+sv)/2 - (ai+si_)/2, 2);   % ROI-average across all 64 pairs
roi_d_val = roi_d_val(~isnan(roi_d_val));

if isempty(roi_d_val)
    error(['s13:no_data -- roi_d_val is empty after processing all subjects. ' ...
        'This means no subject produced usable ImC values -- check the console ' ...
        'output above for "Missing files for..." warnings (wrong file path/key ' ...
        'format) or silent all-NaN segc (channel-label mismatch in fIdx/pIdx, or ' ...
        'no trials matched cueArr/validArr for a given condition). Fix upstream ' ...
        'before re-running -- do not proceed with n=0.']);
end

r0 = local_eval_band(roi_d_val, Delta, r_prior);
fprintf('\nROI-average ImC, d_val, %s: n=%d dz=%.3f t=%.3f p_freq=%.4f BF10=%.3f BF01=%.3f p_TOST=%.4f equivalent=%d\n', ...
    bn, r0{1}, r0{2}, r0{3}, r0{4}, r0{5}, r0{6}, r0{7}, r0{8});
Troi = cell2table({bn, r0{1}, r0{2}, r0{3}, r0{4}, r0{5}, r0{6}, r0{7}, r0{8}}, ...
    'VariableNames', {'band','n','dz','t','p_freq','BF10','BF01','p_TOST','equivalent'});
writetable(Troi, fullfile(cfg.deriv_dir, sprintf('imag_coherence_roi_average_%s.csv', bn)));
fprintf('Saved imag_coherence_roi_average_%s.csv\n', bn);

% ======================= per-pair d_val on ImC, at F2-P4 specifically =======================
d_val_pair = (av+sv)/2 - (ai+si_)/2;   % [subject x pair]
d_target = d_val_pair(:,target_pair); d_target = d_target(~isnan(d_target));
rt = local_eval_band(d_target, Delta, r_prior);
fprintf('F2-P4 ImC, d_val, %s: n=%d dz=%.3f t=%.3f p_freq=%.4f BF10=%.3f BF01=%.3f p_TOST=%.4f equivalent=%d\n', ...
    bn, rt{1}, rt{2}, rt{3}, rt{4}, rt{5}, rt{6}, rt{7}, rt{8});

Tout = table(labs(:), 'VariableNames', {'pair'});
Tout.dz = nan(nF*nT,1); Tout.t = nan(nF*nT,1); Tout.p_freq = nan(nF*nT,1);
Tout.BF10 = nan(nF*nT,1); Tout.BF01 = nan(nF*nT,1); Tout.p_TOST = nan(nF*nT,1); Tout.equivalent = false(nF*nT,1);
for k = 1:nF*nT
    dk = d_val_pair(:,k); dk = dk(~isnan(dk));
    if numel(dk) < 3, continue; end
    rk = local_eval_band(dk, Delta, r_prior);
    Tout.dz(k)=rk{2}; Tout.t(k)=rk{3}; Tout.p_freq(k)=rk{4};
    Tout.BF10(k)=rk{5}; Tout.BF01(k)=rk{6}; Tout.p_TOST(k)=rk{7}; Tout.equivalent(k)=rk{8};
end
writetable(Tout, fullfile(cfg.deriv_dir, sprintf('imag_coherence_valid_main_%s.csv', bn)));
fprintf('\nSaved imag_coherence_valid_main_%s.csv (per-pair BF/TOST on ImC, d_val contrast)\n', bn);

% ======================= volume-conduction diagnostic: phase-difference concentration =======================
Tphase = struct2table(phase_diag);
writetable(Tphase, fullfile(cfg.deriv_dir, sprintf('phase_diff_diagnostic_F2P4_%s.csv', bn)));
fprintf('Saved phase_diff_diagnostic_F2P4_%s.csv\n', bn);
fprintf(['\nInterpretation guide:\n' ...
  '  - If F2-P4 ImC d_val is still significant/BF10>3: the effect likely reflects\n' ...
  '    genuine non-zero-lag connectivity, not volume conduction alone.\n' ...
  '  - If F2-P4 ImC d_val collapses to null (BF01>3, equivalent=true): the PLV\n' ...
  '    effect may be driven substantially by zero-lag (volume conduction/shared\n' ...
  '    reference) coupling and should be reported/discussed with that caveat.\n' ...
  '  - In phase_diff_diagnostic_*.csv, mean_angle_deg near 0 or +-180 with a large\n' ...
  '    resultant_len is itself a volume-conduction signature (phase differences\n' ...
  '    clustering at zero/pi lag); values elsewhere support a true lagged coupling.\n']);

% ======================= local functions =======================
function [set_file, trials_file] = local_subject_files(cfg, key)
    set_file    = fullfile(cfg.deriv_dir, sprintf('%s_clean.set', key));
    trials_file = fullfile(cfg.deriv_dir, sprintf('%s_trials.mat', key));
end

function row = local_eval_band(d, Delta, r_prior)
    n  = numel(d);
    dz = mean(d)/std(d);
    se = std(d)/sqrt(n);
    tobs = mean(d)/se;
    v = n - 1;
    p_freq = 2*(1 - tcdf(abs(tobs), v));

    BF10 = local_jzs_bf10(tobs, n, r_prior);
    BF01 = 1/BF10;

    bound_raw = Delta * std(d);
    t_lo = (mean(d) - (-bound_raw)) / se;
    t_hi = (mean(d) -   bound_raw ) / se;
    p_lo = 1 - tcdf(t_lo, v);
    p_hi = tcdf(t_hi, v);
    p_tost = max(p_lo, p_hi);
    equiv = p_tost < 0.05;

    row = {n, dz, tobs, p_freq, BF10, BF01, p_tost, equiv};
end

function bf10 = local_jzs_bf10(t, n, r)
    % Rouder et al. (2009), one-sample JZS Bayes factor. Copied verbatim
    % from s12_bayes_tost.m -- do not reimplement independently.
    v = n - 1;
    num = (1 + t^2/v)^(-(v+1)/2);
    integrand = @(g) (1+n.*g).^(-0.5) .* ...
        (1 + t^2 ./ ((1+n.*g).*v)).^(-(v+1)/2) .* ...
        (r/sqrt(2*pi)) .* g.^(-1.5) .* exp(-(r^2)./(2.*g));
    denom = integral(integrand, 0, Inf);
    bf01 = num / denom;
    bf10 = 1 / bf01;
end
