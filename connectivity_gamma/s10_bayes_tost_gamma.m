% s10_bayes_tost_gamma.m -- Phase 2: Bayes Factor (JZS) and TOST equivalence
% tests for frontoparietal dwPLI/PLV in the gamma band, Pilot 2 (ds003702).
%
% Mirrors s12_bayes_tost.m (Pilot 1, ds004117) EXACTLY in method:
%   - Delta   = 0.4     (equivalence bound, Cohen's dz)
%   - r_prior = 0.707   (Cauchy prior scale)
%   - jzs_bf10() copied verbatim from s12_bayes_tost.m (manual numerical
%     integration over g, Rouder et al. 2009) -- not from a package.
%   - TOST logic copied verbatim from s12's eval_band(): two one-sided tests
%     against bound_raw = Delta*std(d) on the raw-difference scale.
%
% Design differs from Pilot 1 (single load factor) because ds003702 is a
% 2x2 within-subject design (cue: avatar/stick x validity: valid/invalid).
% Per-subject contrast scores are built BEFORE calling eval_band(), using
% the identical construction already used in s04_stats.m for cue_main /
% valid_main / interaction:
%   d_cue = mean(avatar_valid,avatar_invalid) - mean(stick_valid,stick_invalid)
%   d_val = mean(avatar_valid,stick_valid)     - mean(avatar_invalid,stick_invalid)
%   d_int = (avatar_valid - avatar_invalid)    - (stick_valid - stick_invalid)
%
% Input: conn_results.mat, R table with fields
%   {subject, band, cue, valid, dwpli, plv, frontal_power, parietal_power}
% -- this schema is taken directly from s04_stats.m's local_getval(), which
% pulls R by (subject, band, cue, valid) and reads a single scalar per
% metric (r.(mt)(1)). This resolves the "gamma results filename/variable
% structure" question flagged previously: it is the SAME conn_results.mat /
% R table used by s04, just filtered to band = 'gamma_low' / 'gamma_high'.
% This assumes cfg.bands already contains gamma_low (30-45 Hz) and
% gamma_high (55-80 Hz) entries, since s04 loops generically over
% fieldnames(cfg.bands) with no band restriction. If cfg.bands does NOT yet
% define gamma_low/gamma_high, conn_results.mat needs to be regenerated
% (s03) with those bands added before this script will find any rows.
%
% Output:
%   1. bayes_tost_gamma.csv -- 12 rows (3 contrasts x 2 metrics x 2 subbands)
%      columns: contrast, metric, band, n, dz, t, p_freq, BF10, BF01,
%               p_TOST, equivalent_at_0p4
%   2. Console table, printed in the same style as s12.
%   3. Per-pair BF01 heatmap -- NOT IMPLEMENTED HERE. See note at bottom.

clear; clc;
cfg = config();
Delta   = 0.4;    % equivalence bound, Cohen's dz -- copied from s12, unchanged
r_prior = 0.707;  % Cauchy prior scale -- copied from s12, unchanged

load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');  % same file/variable as s04

metrics = {'dwpli','plv'};
subbands = {'gamma_low','gamma_high'};
subs = unique(R.subject);

% helper: pull one scalar per subject/band/metric/cue/valid, or NaN.
% Identical in form to local_getval() in s04_stats.m.
getval = @(subj,bn,mt,cueLab,valLab) local_getval(R,subj,bn,mt,cueLab,valLab);

rows = {};
fprintf('\n contrast     metric    band          n     dz      t       p_freq    BF10      BF01     p_TOST   equiv\n');
fprintf('----------------------------------------------------------------------------------------------------------\n');

for s = 1:numel(subbands)
    bn = subbands{s};
    for m = 1:numel(metrics)
        mt = metrics{m};

        av = nan(numel(subs),1); ai = nan(numel(subs),1);
        sv = nan(numel(subs),1); si = nan(numel(subs),1);
        for i = 1:numel(subs)
            av(i) = getval(subs{i}, bn, mt, 'avatar', 'valid');
            ai(i) = getval(subs{i}, bn, mt, 'avatar', 'invalid');
            sv(i) = getval(subs{i}, bn, mt, 'stick',  'valid');
            si(i) = getval(subs{i}, bn, mt, 'stick',  'invalid');
        end

        % per-subject contrast scores, built exactly as in s04_stats.m
        contrasts = struct();
        contrasts.d_cue = (av+ai)/2 - (sv+si)/2;
        contrasts.d_val = (av+sv)/2 - (ai+si)/2;
        contrasts.d_int = (av-ai)   - (sv-si);

        cnames = fieldnames(contrasts);
        for c = 1:numel(cnames)
            cn = cnames{c};
            d = contrasts.(cn);
            d = d(~isnan(d));

            r0 = eval_band(bn, mt, d, Delta, r_prior);  % {band,metric,n,dz,t,p_freq,BF10,BF01,p_TOST,equiv}
            row = [{cn}, r0(2), r0(1), r0(3:end)];       % reorder to contrast,metric,band,...
            rows(end+1,:) = row; %#ok<SAGROW>

            fprintf('%-10s   %-8s  %-12s %3d  %6.3f  %6.2f   %.4f    %6.3f   %6.3f   %.4f   %d\n', ...
                row{1}, row{2}, row{3}, row{4}, row{5}, row{6}, row{7}, row{8}, row{9}, row{10}, row{11});
        end
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'contrast','metric','band','n','dz','t','p_freq','BF10','BF01','p_TOST','equivalent_at_0p4'});
disp(T);
writetable(T, fullfile(cfg.deriv_dir,'bayes_tost_gamma.csv'));
fprintf('\nSaved bayes_tost_gamma.csv -- gamma-band Bayes/TOST table (3 contrasts x 2 metrics x 2 subbands)\n');

% ======================= local functions (copied verbatim from s12) =======================
function v = local_getval(R, subj, bn, mt, cueLab, valLab)
    r = R(strcmp(R.subject,subj) & strcmp(R.band,bn) & ...
          strcmp(R.cue,cueLab) & strcmp(R.valid,valLab), :);
    if isempty(r), v = NaN; else, v = r.(mt)(1); end
end

function row = eval_band(bn, mt, d, Delta, r_prior)
    n  = numel(d);
    dz = mean(d)/std(d);
    se = std(d)/sqrt(n);
    tobs = mean(d)/se;
    v = n - 1;
    p_freq = 2*(1 - tcdf(abs(tobs), v));

    BF10 = jzs_bf10(tobs, n, r_prior);
    BF01 = 1/BF10;

    bound_raw = Delta * std(d);
    t_lo = (mean(d) - (-bound_raw)) / se;
    t_hi = (mean(d) -   bound_raw ) / se;
    p_lo = 1 - tcdf(t_lo, v);
    p_hi = tcdf(t_hi, v);
    p_tost = max(p_lo, p_hi);
    equiv = p_tost < 0.05;

    row = {bn, mt, n, dz, tobs, p_freq, BF10, BF01, p_tost, equiv};
end

function bf10 = jzs_bf10(t, n, r)
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

% ======================= NOT IMPLEMENTED: per-pair heatmap =======================
% Item 2 of the request (128-cell BF01 heatmap: 64 frontal x parietal pairs
% x 2 subbands, for each of the 3 contrasts) requires per-pair-level
% dwPLI/PLV values per subject/condition/subband. The R table used above
% (and in s04_stats.m) only carries ROI-AVERAGED values -- one scalar per
% subject/band/cue/valid per metric (r.(mt)(1)) -- not per-channel-pair
% values. No upstream script or .mat file confirming a per-pair-level gamma
% output (filename, variable name, or table schema with a pair identifier)
% has been reviewed yet. Per the "schema verification before coding"
% protocol, this is flagged rather than guessed: before the heatmap can be
% written, confirm (a) which script produces per-pair connectivity for the
% gamma bands, (b) its output filename, and (c) whether it stores one row
% per pair or a pair x subband x subject array. Please point me to that
% script/file (or confirm it doesn't exist yet and needs to be added
% upstream), and I'll add the heatmap as a follow-up (e.g. s11_gamma_heatmap.m).
