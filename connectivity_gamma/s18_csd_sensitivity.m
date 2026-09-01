% s18_csd_sensitivity.m -- Sensitivity branch: CSD (Surface Laplacian,
% Perrin et al. 1989) re-reference, as a differential test between two
% mechanisms that dwPLI/Imaginary Coherency alone cannot distinguish:
%   (a) genuine volume conduction (shared cortical/skull spread), vs.
%   (b) an average-reference artifact specific to THIS pipeline's
%       pop_reref(EEG, []) step (s01_preprocess_New.m).
% A reviewer raised this as a substantive limitation of the PLV/dwPLI/ImC
% dissociation reported for valid_main / gamma_high (S10/S13).
%
% WHY THIS BRANCHES FROM THE ALREADY-SAVED <key>_clean.set (NOT FROM RAW):
% s01_preprocess_New.m applies average reference as its LAST step, after
% pop_interp. No intermediate (post-interp, pre-reref) file was saved, so
% branching there would require re-running filtering/clean_rawdata/ICA/
% ICLabel for all 47 subjects -- expensive and unnecessary. This is safe to
% skip because the surface Laplacian is mathematically REFERENCE-FREE (a
% spatial second-derivative operator that cancels any common additive
% reference term by construction; Perrin et al. 1989; Kayser & Tenke, 2006).
% Applying ft_scalpcurrentdensity to the average-referenced clean.set is
% therefore numerically equivalent to applying it before pop_reref. This
% equivalence is the justification for using the existing final files
% rather than re-deriving a pre-reref intermediate.
%
% WHY FieldTrip ft_scalpcurrentdensity('spline') AND NOT THE CSD TOOLBOX
% (Kayser & Tenke): FieldTrip is already an installed, path-configured
% dependency in this pipeline (s08/s09 alpha validation), and
% eeglab2fieldtrip(EEG,'raw','none') already produces a valid 3-D `elec`
% structure (confirmed via s08_patch_add_elec.m, which extracts exactly
% this field) -- no new montage/coordinate file to source or verify.
% FieldTrip's 'spline' method implements the same Perrin et al. (1989)
% algorithm the CSD Toolbox uses, so there is no methodological advantage
% to switching, only the added risk of an unverified montage-specific
% .csd file (the same class of risk that caused the earlier
% neighbours.mat/layout mixup in s06c).
%
% CHANNEL-SET CONSISTENCY: EEG.urchanlocs in s01 is captured BEFORE the
% per-subject manual bad-channel removal (e.g. sub-48's P1), so the final
% pop_interp step restores it -- every subject's clean.set has the same
% full 61-channel montage, satisfying ft_scalpcurrentdensity's requirement
% of a consistent channel set.
%
% SCOPE (per request): valid_main / gamma_high only, at two levels:
%   (1) ROI-average (all 64 frontal-parietal pairs, mean dwPLI/PLV/ImC)
%   (2) F2-P4 pair alone
% Metrics: dwPLI, PLV, Imaginary Coherency (ImC) -- same formulas as
% s03_connectivity.m / s13_imaginary_coherence.m, applied to CSD-
% transformed data instead of average-referenced data. NOT using
% ft_connectivityanalysis/ft_freqstatistics, to stay equivalent to the
% pipeline's hand-written Hilbert-based connectivity code elsewhere.
%
% Output:
%   csd_sensitivity_summary.csv    -- group-level BF/TOST/frequentist
%                                      table, one row per level x metric
%   csd_sensitivity_persubject.csv -- per-subject raw d_val contrast
%                                      values, for transparency/reuse
%                                      (BF/t/p are group-level only and
%                                      are NOT duplicated per subject here)

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;

M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

bn = 'gamma_high';
band_range = cfg.bands.(bn);
srate = cfg.resample_hz;
win = round(cfg.maint_window * srate);
wlen = win(2) - win(1) + 1;

nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_parietal);
labs = cell(1,nF*nT);
for i=1:nF, for j=1:nT, labs{(i-1)*nT+j} = sprintf('%s-%s', cfg.roi_frontal{i}, cfg.roi_parietal{j}); end, end
target_pair = find(strcmp(labs,'F2-P4'));

cueLabs = {'avatar','stick'}; cueVals = [cfg.ti.cue_avatar, cfg.ti.cue_stick];
valLabs = {'valid','invalid'}; valVals = [cfg.ti.val_valid, cfg.ti.val_invalid];

% per-subject, per-condition, per-pair storage for all 3 metrics
nSub = height(manifest);
DWPLI = nan(nSub, 4, nF*nT); PLV = nan(nSub, 4, nF*nT); IMC = nan(nSub, 4, nF*nT);
subj_list = cell(nSub,1);
n_ok = 0;

for su = 1:nSub
    subj = char(manifest.subject(su)); key = char(manifest.key(su));
    subj_list{su} = subj;

    trials_file = fullfile(cfg.deriv_dir, [key '_trials.mat']);
    set_file = fullfile(cfg.deriv_dir, [key '_clean.set']);
    if ~isfile(trials_file) || ~isfile(set_file)
        warning('%s: missing trials/clean.set, skipping.', subj); continue;
    end
    T = load(trials_file); trials = T.trials;
    good_trials = trials([trials.good]);
    if isempty(good_trials), warning('%s: no good trials, skipping.', subj); continue; end

    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);

    % ---- CSD re-reference (branch point; see header note on equivalence) ----
    ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
    cfg_csd = []; cfg_csd.method = 'spline';
    ft_csd = ft_scalpcurrentdensity(cfg_csd, ft_raw);
    csd_data = ft_csd.trial{1};             % [nchan x ntime], CSD-transformed continuous data
    csd_labels = ft_csd.label;

    % channel index lookup against CSD output labels (case-sensitive match,
    % same channel set/order expected as EEG.chanlocs since CSD doesn't
    % reorder or drop channels)
    fr = local_chanidx(csd_labels, cfg.roi_frontal);
    pa = local_chanidx(csd_labels, cfg.roi_parietal);
    if numel(fr) < numel(cfg.roi_frontal) || numel(pa) < numel(cfg.roi_parietal)
        warning('%s: missing ROI channel(s) in CSD output, skipping.', subj); continue;
    end

    Xf = util_bandpass(double(csd_data), srate, band_range, cfg.filt_order_factor);
    Z  = hilbert(Xf.').';

    cond_ok = true(1,4);
    for c = 1:2
        for v = 1:2
            cidx = (c-1)*2 + v;
            sel = good_trials([good_trials.cue]==cueVals(c) & [good_trials.valid]==valVals(v));
            ntr = numel(sel);
            if ntr < 5, cond_ok(cidx) = false; continue; end

            seg = nan(numel(csd_labels), wlen, ntr); kept = 0;
            for t = 1:ntr
                a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                if a < 1 || z > size(Z,2), continue; end
                kept = kept+1; seg(:,:,kept) = Z(:, a:z);
            end
            seg = seg(:,:,1:kept);
            if kept < 5, cond_ok(cidx) = false; continue; end

            p = 0;
            for i = 1:numel(fr)
                zi = reshape(seg(fr(i),:,:), wlen, kept);
                for j = 1:numel(pa)
                    p = p+1;
                    zj = reshape(seg(pa(j),:,:), wlen, kept);
                    X = zi.*conj(zj); im = imag(X);
                    s_im=sum(im,2); s_abs=sum(abs(im),2); s_sq=sum(im.^2,2);
                    dwpli_t = (s_im.^2 - s_sq) ./ (s_abs.^2 - s_sq);
                    plv_t = abs(mean(X./abs(X), 2));
                    DWPLI(su,cidx,p) = mean(dwpli_t,'omitnan');
                    PLV(su,cidx,p) = mean(plv_t,'omitnan');

                    Sxy = mean(mean(X)); Sxx = mean(mean(abs(zi).^2)); Syy = mean(mean(abs(zj).^2));
                    IMC(su,cidx,p) = imag(Sxy/sqrt(Sxx*Syy));
                end
            end
        end
    end
    n_ok = n_ok + 1;
    fprintf('%s: CSD-based connectivity computed (%d/4 conditions usable).\n', subj, sum(cond_ok));
end
fprintf('\nSubjects processed: %d/%d.\n', n_ok, nSub);

% ======================= build d_val contrast, both levels, all 3 metrics =======================
% condition index: 1=avatar_valid, 2=avatar_invalid, 3=stick_valid, 4=stick_invalid
Delta = 0.4; r_prior = 0.707;   % copied verbatim from s10/s12/s13
metrics = struct('dwPLI', DWPLI, 'PLV', PLV, 'ImC', IMC);
mnames = fieldnames(metrics);
levels = {'roi_avg','F2P4'};

summary_rows = {}; persubj_rows = {};
for mi = 1:numel(mnames)
    Mval = metrics.(mnames{mi});
    av = Mval(:,1,:); ai = Mval(:,2,:); sv = Mval(:,3,:); iv = Mval(:,4,:);
    d_val_pair = squeeze((av+sv)/2 - (ai+iv)/2);   % [subject x pair]

    for li = 1:numel(levels)
        if strcmp(levels{li},'roi_avg')
            d = mean(d_val_pair, 2, 'omitnan');
        else
            d = d_val_pair(:, target_pair);
        end
        for si = 1:nSub
            if ~isnan(d(si))
                persubj_rows(end+1,:) = {subj_list{si}, levels{li}, mnames{mi}, d(si)}; %#ok<SAGROW>
            end
        end
        d = d(~isnan(d)); n = numel(d);
        if n < 3
            warning('%s / %s: n=%d, too few subjects for group stats.', levels{li}, mnames{mi}, n);
            continue;
        end
        dz = mean(d)/std(d); se = std(d)/sqrt(n); tobs = mean(d)/se; v = n-1;
        p_freq = 2*(1-tcdf(abs(tobs),v));
        BF10 = local_jzs_bf10(tobs, n, r_prior); BF01 = 1/BF10;
        bound_raw = Delta*std(d);
        t_lo=(mean(d)-(-bound_raw))/se; t_hi=(mean(d)-bound_raw)/se;
        p_tost = max(1-tcdf(t_lo,v), tcdf(t_hi,v));
        summary_rows(end+1,:) = {levels{li}, mnames{mi}, n, dz, tobs, p_freq, BF10, BF01, p_tost, p_tost<0.05}; %#ok<SAGROW>
    end
end

Tsum = cell2table(summary_rows, 'VariableNames', ...
    {'level','metric','n','dz','t_freq','p_freq','BF10','BF01','p_TOST','equivalent_at_0p4'});
disp(Tsum);
writetable(Tsum, fullfile(cfg.deriv_dir, 'csd_sensitivity_summary.csv'));

Tps = cell2table(persubj_rows, 'VariableNames', {'subject','level','metric','d_val'});
writetable(Tps, fullfile(cfg.deriv_dir, 'csd_sensitivity_persubject.csv'));

fprintf('\nSaved csd_sensitivity_summary.csv and csd_sensitivity_persubject.csv\n');
fprintf(['\nInterpretation guide:\n' ...
  '  - Compare PLV row (roi_avg and F2P4) here to the average-reference results\n' ...
  '    (S10: BF10=45.16 roi_avg; F2-P4 max-stat t=-3.70).\n' ...
  '  - If the PLV effect COLLAPSES under CSD (BF01>3, small |t_freq|): the original\n' ...
  '    average-reference PLV effect was likely a reference artifact, not volume\n' ...
  '    conduction alike -- CSD removes both reference artifacts AND genuine\n' ...
  '    volume-conducted zero-lag spread, so it alone cannot fully separate the\n' ...
  '    two, but it rules out reference choice as a confound either way.\n' ...
  '  - dwPLI/ImC are expected to remain null under CSD too (both already immune to\n' ...
  '    zero-lag under either reference) -- this is a consistency check, not a new test.\n']);

% ======================= local functions =======================
function idx = local_chanidx(labels, want)
    idx = [];
    for i = 1:numel(want)
        j = find(strcmpi(labels, want{i}), 1);
        if isempty(j)
            warning('channel "%s" not found in CSD output labels', want{i});
        else
            idx(end+1) = j; %#ok<AGROW>
        end
    end
end

function bf10 = local_jzs_bf10(t, n, r)
    % Rouder et al. (2009), one-sample JZS Bayes factor. Copied verbatim
    % from s12_bayes_tost.m -- do not reimplement independently.
    v = n-1; num = (1+t^2/v)^(-(v+1)/2);
    integrand = @(g) (1+n.*g).^(-0.5).*(1+t^2./((1+n.*g).*v)).^(-(v+1)/2).*(r/sqrt(2*pi)).*g.^(-1.5).*exp(-(r^2)./(2.*g));
    denom = integral(integrand,0,Inf); bf10 = 1/(num/denom);
end
