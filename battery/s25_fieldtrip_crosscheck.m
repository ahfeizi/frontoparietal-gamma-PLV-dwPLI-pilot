% s25_fieldtrip_crosscheck.m -- Estimator VALIDITY check (not reproducibility).
% Our PLV, dwPLI and ImC (util_conn_pairs.m, the s21_battery.m path) are
% compared with independent implementations computed from the SAME analytic
% signal:
%   subject sub-01 | band gamma_high | band-pass + Hilbert | average reference
%   all four cue x validity cells (valid and invalid under each cue) | all 64
%   frontal x parietal pairs
%
% Analytic signal: as s21_battery.m -- data loading, good trials, ROI rows,
% util_bandpass + hilbert on the continuous ROI data, maintenance window
% locked to trialinfo Maint, out-of-bounds trials dropped, cells with < 5
% segments skipped. The local functions are copied from s21 (local functions
% cannot be shared between scripts).
%
% Independent estimates, per pair (zi, zj = [samples x trials] analytic signal):
%   dwPLI     FieldTrip ft_connectivity_wpli, debias = true. Input is the
%             analytic signal as Fourier coefficients [trials x chan x samples]
%             with isunivariate = true, so FieldTrip forms the cross-spectra
%             itself (as ft_connectivityanalysis does for fourierspctrm input).
%             Output = dwPLI across trials at each sample; then mean over
%             samples (omitnan) -- the aggregation of util_conn_pairs.
%   ImC       FieldTrip ft_checkdata(..., 'cmbstyle', 'fullfast') converts the
%             coefficients, every sample x trial one observation, into the
%             cross-spectral matrix averaged jointly over samples x trials
%             (the conversion ft_connectivityanalysis uses); then
%             ft_connectivity_corr(..., 'complex', 'imag', 'pownorm', true)
%             = imag(Sxy / sqrt(Sxx*Syy)). Same aggregation as ours (joint).
%   PLV       directly from dphi = angle(zi) - angle(zj):
%             PLV       mean over samples of |mean over trials exp(1i*dphi)|.
%                       This is the battery's PLV aggregation (util_conn_pairs,
%                       verbatim from s03): across trials per sample, then
%                       averaged over samples.
%             PLV_joint |mean(exp(1i*dphi))| over samples x trials. NOT the
%                       battery's PLV measure; it is resultant_len of
%                       battery_phasediff.csv, so it is compared with
%                       abs(mvec) from util_conn_pairs. Note that util_conn_pairs
%                       computes mvec with the same formula, so this line is a
%                       consistency check, not an independent one.
%
% Tolerance: 1e-10 absolute for every measure. FieldTrip does not force a
% looser one: all estimators are closed-form in double precision and differ
% from ours only in summation order and BLAS vs element-wise products
% (expected max |diff| ~ 1e-15 .. 1e-13).
%
% Also checks that our values equal the sub-01 rows of battery_persubject.csv
% in cfg.battery.check_dir (confirms the analytic signal is the one s21 used;
% reported on its own line, outside the estimator OVERALL).
%
% PATH: the only FieldTrip on the path is cfg.fieldtrip_dir, which must be the
% full FieldTrip copy in the Documents folder. After `eeglab nogui`, EEGLAB's
% bundled Fieldtrip-lite / Fileio plugin folders are removed from the path,
% then ft_defaults is called, and every FieldTrip function used is checked to
% resolve inside cfg.fieldtrip_dir only (error otherwise).
%
% Writes nothing. Run with MATLAB's current folder at the repository root
% (config.m paths are relative to it).

clear; clc;
repo = fileparts(fileparts(mfilename('fullpath')));   % repository root (this script is in battery/)
addpath(repo, fullfile(repo, 'battery'));             % config.m, util_bandpass.m, util_chanidx.m; util_conn_pairs.m
cfg = config();

subj_target = 'sub-01';
bn   = 'gamma_high';
band = cfg.battery.bands.(bn);
tol  = 1e-10;

% ---------------- path: EEGLAB without its FieldTrip copy, then FieldTrip ----------------
eeglab_dir = local_abspath(cfg.eeglab_dir);
ft_dir     = local_abspath(cfg.fieldtrip_dir);
if ~isfile(fullfile(ft_dir, 'ft_defaults.m'))
    error('cfg.fieldtrip_dir resolves to %s, which has no ft_defaults.m.', ft_dir);
end
if local_under(ft_dir, eeglab_dir)
    error('cfg.fieldtrip_dir (%s) is inside the EEGLAB tree; point it at the FieldTrip copy in Documents.', ft_dir);
end
if ~contains(lower([ft_dir filesep]), lower([filesep 'Documents' filesep]))
    error('cfg.fieldtrip_dir resolves to %s, not to the FieldTrip copy in the Documents folder.', ft_dir);
end

addpath(eeglab_dir); eeglab nogui;
removed = local_strip_ft(eeglab_dir);
fprintf('Removed %d EEGLAB FieldTrip-lite/Fileio folder(s) from the path.\n', numel(removed));
if ~isempty(removed), fprintf('  %s\n', removed{:}); end
% drop any FieldTrip state set up by the lite copy before initialising ours
clear global ft_default
clear('ft_defaults','ft_hastoolbox','ft_platform_supports','ft_getopt','ft_warning','ft_notice', ...
      'ft_info','ft_debug','ft_error','ft_progress','ft_checkdata','ft_datatype','ft_datatype_freq', ...
      'ft_connectivity_wpli','ft_connectivity_corr');
addpath(ft_dir); ft_defaults;
addpath(fullfile(ft_dir, 'connectivity'));            % ft_connectivity_*; ft_defaults does not add it
local_check_path(ft_dir, eeglab_dir, ...
    {'ft_defaults','ft_checkdata','ft_datatype_freq','ft_getopt','ft_progress', ...
     'ft_connectivity_wpli','ft_connectivity_corr'}, {'nansum','hilbert','fir1','filtfilt'});
fprintf('FieldTrip in use: %s\n', fileparts(which('ft_defaults')));
% the analytic signal must use the same filter/Hilbert code as s21 (MATLAB Signal
% Processing Toolbox, not FieldTrip's external/signal drop-ins)
for f = {'hilbert','fir1','filtfilt'}
    w = which(f{1});
    fprintf('  %-8s %s\n', f{1}, w);
    if ~local_under(w, matlabroot)
        warning('%s does not come from the MATLAB installation; the analytic signal may differ from s21.', f{1});
    end
end

% ---------------- data: as s21 ----------------
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;
su = find(strcmp(cellstr(string(manifest.subject)), subj_target), 1);
if isempty(su), error('%s not found in manifest.mat', subj_target); end
key = char(manifest.key(su));

srate = cfg.resample_hz;
cueLabs = {'avatar','stick'}; cueVals = [cfg.ti.cue_avatar, cfg.ti.cue_stick];
valLabs = {'valid','invalid'}; valVals = [cfg.ti.val_valid, cfg.ti.val_invalid];
conds = {}; condCue = []; condVal = [];
for c = 1:2
    for v = 1:2
        conds{end+1} = sprintf('%s_%s', cueLabs{c}, valLabs{v}); %#ok<SAGROW>
        condCue(end+1) = cueVals(c); condVal(end+1) = valVals(v); %#ok<SAGROW>
    end
end
nC = numel(conds);

EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
T = load(fullfile(cfg.deriv_dir, [key '_trials.mat'])); trials = T.trials;
good_trials = trials([trials.good]);
if isempty(good_trials), error('%s: no good trials.', subj_target); end

have = {EEG.chanlocs.labels};
labF = local_present(have, cfg.roi_frontal);
labP = local_present(have, cfg.roi_parietal);
if numel(labF) < numel(cfg.roi_frontal) || numel(labP) < numel(cfg.roi_parietal)
    warning('%s: only %d/%d frontal and %d/%d parietal ROI channels found -- check labels.', ...
        subj_target, numel(labF), numel(cfg.roi_frontal), numel(labP), numel(cfg.roi_parietal));
end
roiLabs = [labF labP];
iAvg = util_chanidx(EEG, roiLabs);

runs = struct('srate', srate, 'sel', {cell(1,nC)}, 'average', double(EEG.data(iAvg,:)));
for c = 1:nC
    runs.sel{c} = good_trials([good_trials.cue] == condCue(c) & [good_trials.valid] == condVal(c));
end
clear EEG

nF = numel(labF); nP = numel(labP); np = nF*nP;
fr = 1:nF; pa = nF + (1:nP);            % indices into the ROI-row segments
pairLabs = cell(1,np);
for i = 1:nF, for j = 1:nP, pairLabs{(i-1)*nP+j} = sprintf('%s-%s', labF{i}, labP{j}); end, end

seg = local_pool(runs, 'average', cfg.maint_window, ...
    @(X, sr) local_hilbert(X, sr, band, cfg.filt_order_factor));

% ---------------- ours vs independent ----------------
meas  = {'PLV','PLV_joint','dwPLI','ImC'};
desc  = {'trials per sample, mean over samples (battery PLV)', ...
         'joint samples x trials vs abs(mvec) (resultant_len)', ...
         'FieldTrip ft_connectivity_wpli, debias', ...
         'FieldTrip fullfast CSD + ft_connectivity_corr imag'};
ours  = nan(nC, np, numel(meas));
indep = nan(nC, np, numel(meas));
for c = 1:nC
    kept = local_ntr(seg{c});
    if kept < 5
        fprintf('  %s %s: only %d in-bounds segments, skipping cell (as s21)\n', subj_target, conds{c}, kept);
        continue;
    end
    fprintf('  %s: %d trials x %d samples\n', conds{c}, kept, size(seg{c},2));
    [plv, dwpli, imc, mvec] = util_conn_pairs(seg{c}, fr, pa);
    ours(c,:,:) = cat(3, plv, abs(mvec), dwpli, imc);

    [plv_d, plv_joint_d] = local_plv_direct(seg{c}, fr, pa);
    indep(c,:,:) = cat(3, plv_d, plv_joint_d, ...
        local_ft_dwpli(seg{c}, fr, pa), ...
        local_ft_imc(seg{c}, fr, pa, roiLabs, mean(band)));
end
done_c = find(~all(isnan(ours(:,:,1)), 2)).';
if isempty(done_c), error('%s: no cell with >= 5 segments.', subj_target); end

fprintf('\n=== %s | %s | hilbert | average | %d cell(s) x %d pairs | tol %.0e ===\n', ...
    subj_target, bn, numel(done_c), np, tol);
nfail = 0;
for m = 1:numel(meas)
    a = ours(done_c,:,m); b = indep(done_c,:,m);
    d = abs(a - b);
    d(isnan(a) & isnan(b)) = 0;
    d(isnan(a) ~= isnan(b)) = Inf;          % NaN on one side only is a failure
    [dmax, k] = max(d(:));
    [ci, p] = ind2sub(size(d), k);
    ok = dmax <= tol;
    if ok, tag = 'PASS'; else, tag = 'FAIL'; nfail = nfail + 1; end
    fprintf('%s  %-9s %-52s max|ours-indep| = %.3e (%s, %s); ours in [%.4g, %.4g]\n', ...
        tag, meas{m}, desc{m}, dmax, conds{done_c(ci)}, pairLabs{p}, min(a(:)), max(a(:)));
end
if nfail == 0
    fprintf('\nOVERALL (estimators): PASS\n');
else
    fprintf('\nOVERALL (estimators): FAIL (%d measure(s) failed)\n', nfail);
end

% ---------------- provenance: our values == s21 output for sub-01 ----------------
f_b = fullfile(cfg.battery.check_dir, 'battery_persubject.csv');
if ~isfile(f_b)
    fprintf('SKIP  s21 provenance: %s not found\n', f_b);
else
    L = splitlines(fileread(f_b));
    L = L(startsWith(L, [subj_target ',']));
    C = cell(0, 8);
    if numel(L) > 1
        C = split(L, ',');                  % subject,condition,band,measure,decomposition,reference,pair,value
    end
    C = C(strcmp(C(:,3), bn) & strcmp(C(:,5), 'hilbert') & strcmp(C(:,6), 'average'), :);
    keyCsv = strcat(C(:,2), '|', C(:,4), '|', C(:,7));
    valCsv = str2double(C(:,8));
    bm = {'PLV','dwPLI','ImC'}; om = [1 3 4];   % battery measure -> column of ours
    dmax = 0; nmiss = 0;
    for c = done_c
        for m = 1:numel(bm)
            keyOurs = strcat(conds{c}, '|', bm{m}, '|', pairLabs(:));
            [tf, loc] = ismember(keyOurs, keyCsv);
            nmiss = nmiss + sum(~tf);
            if any(tf)
                dmax = max(dmax, max(abs(reshape(ours(c,tf,om(m)), [], 1) - valCsv(loc(tf)))));
            end
        end
    end
    if nmiss == 0 && dmax <= tol
        fprintf('PASS  s21 provenance: ours == %s rows of %s, max|diff| = %.3e\n', subj_target, f_b, dmax);
    else
        fprintf('FAIL  s21 provenance vs %s: %d value(s) missing, max|diff| = %.3e\n', f_b, nmiss, dmax);
    end
end

% ======================= local functions =======================
function [plv, plv_joint] = local_plv_direct(seg, fr, pa)
    % PLV straight from the phase difference of the analytic signal
    wlen = size(seg,2); kept = size(seg,3);
    ph = angle(seg);
    np = numel(fr) * numel(pa); plv = zeros(1,np); plv_joint = zeros(1,np); p = 0;
    for i = 1:numel(fr)
        for j = 1:numel(pa)
            p = p + 1;
            dphi = reshape(ph(fr(i),:,:) - ph(pa(j),:,:), wlen, kept);   % [samples x trials]
            e = exp(1i*dphi);
            plv(p)       = mean(abs(mean(e, 2)), 'omitnan');            % across trials, then over samples
            plv_joint(p) = abs(mean(e(:)));                             % jointly over samples x trials
        end
    end
end

function dw = local_ft_dwpli(seg, fr, pa)
    % FieldTrip debiased wPLI across trials at each sample, then mean over samples
    F = permute(seg, [3 1 2]);                                         % [trials x chan x samples]
    kept = size(F,1);
    w = ft_connectivity_wpli(F, 'debias', true, 'dojack', false, ...
        'isunivariate', true, 'cumtapcnt', ones(kept,1));              % [chan x chan x samples]
    if any(imag(w(:)) ~= 0), error('ft_connectivity_wpli returned a non-real dwPLI.'); end
    w = real(w);
    np = numel(fr) * numel(pa); dw = zeros(1,np); p = 0;
    for i = 1:numel(fr)
        for j = 1:numel(pa)
            p = p + 1;
            dw(p) = mean(reshape(w(fr(i),pa(j),:), [], 1), 'omitnan');
        end
    end
end

function imc = local_ft_imc(seg, fr, pa, labels, fcentre)
    % FieldTrip cross-spectrum averaged jointly over samples x trials, then imag coherency
    [nR, wlen, kept] = size(seg);
    F = reshape(permute(seg, [2 3 1]), wlen*kept, nR);                 % [samples*trials x chan]
    % fullfast divides by the number of non-zero coefficients of channel 1
    if any(F(:,1) == 0), error('zero analytic coefficient in channel 1; fullfast divisor would differ.'); end
    freq = [];
    freq.label         = labels(:);
    freq.freq          = fcentre;                                      % label only, not used in the estimate
    freq.dimord        = 'rpttap_chan_freq';
    freq.fourierspctrm = F;
    freq.cumtapcnt     = repmat(wlen, kept, 1);
    freq = ft_checkdata(freq, 'cmbstyle', 'fullfast');                 % crsspctrm = F.'*conj(F) / n
    crs  = reshape(freq.crsspctrm, [1 nR nR]);
    c = ft_connectivity_corr(crs, 'dimord', 'rpt_chan_chan', 'complex', 'imag', ...
        'pownorm', true, 'hasjack', false);
    c = reshape(c, nR, nR);
    np = numel(fr) * numel(pa); imc = zeros(1,np); p = 0;
    for i = 1:numel(fr)
        for j = 1:numel(pa)
            p = p + 1;
            imc(p) = c(fr(i), pa(j));
        end
    end
end

function Z = local_hilbert(X, srate, band, order_factor)
    % as s21_battery.m / s03_connectivity.m
    Xf = util_bandpass(X, srate, band, order_factor);
    Z  = hilbert(Xf.').';
end

function seg = local_pool(runs, rf, maint_window, transform)
    % as s21_battery.m
    seg = cell(1, numel(runs(1).sel));
    for k = 1:numel(runs)
        srate = runs(k).srate; win = round(maint_window*srate); wlen = win(2)-win(1)+1;
        Z = transform(runs(k).(rf), srate);
        for c = 1:numel(seg)
            seg{c} = cat(3, seg{c}, local_segments(Z, runs(k).sel{c}, win, wlen));
        end
    end
end

function seg = local_segments(Z, sel, win, wlen)
    % as s21_battery.m
    seg = complex(zeros(size(Z,1), wlen, numel(sel))); kept = 0;
    for t = 1:numel(sel)
        a = sel(t).lock_sample + win(1); z = a + wlen - 1;
        if a < 1 || z > size(Z,2), continue; end
        kept = kept + 1;
        seg(:,:,kept) = Z(:, a:z);
    end
    seg = seg(:,:,1:kept);
end

function n = local_ntr(seg)
    if isempty(seg), n = 0; else, n = size(seg,3); end
end

function labs = local_present(have, want)
    % ROI labels found in the montage, in config order (what util_chanidx returns)
    labs = want(cellfun(@(l) any(strcmpi(have, l)), want));
end

function p = local_abspath(p)
    % absolute, canonical form of a (possibly relative) folder
    if ~(startsWith(p, {'\', '/'}) || (numel(p) >= 2 && p(2) == ':'))
        p = fullfile(pwd, p);
    end
    p = char(java.io.File(p).getCanonicalPath());
end

function tf = local_under(p, root)
    p    = [regexprep(p,    '[\\/]+$', '') filesep];
    root = [regexprep(root, '[\\/]+$', '') filesep];
    tf = strncmpi(strrep(p, '/', filesep), strrep(root, '/', filesep), numel(root));
end

function removed = local_strip_ft(eeglab_dir)
    % remove EEGLAB's bundled FieldTrip (plugins/Fieldtrip-lite*, plugins/Fileio*)
    plug = fullfile(eeglab_dir, 'plugins');
    P = strsplit(path, pathsep);
    bad = false(size(P));
    for k = 1:numel(P)
        if ~local_under(P{k}, plug), continue; end
        top = strtok(P{k}(numel(plug)+2:end), '\/');
        bad(k) = ~isempty(regexpi(top, '^(fieldtrip|fileio)', 'once'));
    end
    removed = P(bad);
    if ~isempty(removed), rmpath(removed{:}); end
end

function local_check_path(ft_dir, eeglab_dir, ft_fns, other_fns)
    % 1) no FieldTrip-like path entry outside ft_dir
    P = strsplit(path, pathsep);
    P = P(~cellfun(@isempty, P));
    % (MATLAB's own toolboxes have unrelated fileIO folders, so matlabroot is exempt)
    stray = P(~cellfun(@(p) local_under(p, ft_dir) || local_under(p, matlabroot), P) & ...
              ~cellfun(@isempty, regexpi(P, 'fieldtrip|[\\/]fileio', 'once')));
    if ~isempty(stray)
        error('FieldTrip folders outside %s are on the path:\n  %s', ft_dir, strjoin(stray, '\n  '));
    end
    % 2) every copy of every FieldTrip function used resolves inside ft_dir
    for k = 1:numel(ft_fns)
        hits = which(ft_fns{k}, '-all');
        if isempty(hits), error('%s not found on the path.', ft_fns{k}); end
        out = hits(~cellfun(@(h) local_under(h, ft_dir), hits));
        if ~isempty(out)
            error('%s resolves outside %s:\n  %s', ft_fns{k}, ft_dir, strjoin(out, '\n  '));
        end
    end
    % 3) helpers that may come from MATLAB toolboxes must not resolve to EEGLAB
    for k = 1:numel(other_fns)
        h = which(other_fns{k});
        if isempty(h), error('%s not found on the path.', other_fns{k}); end
        if local_under(h, eeglab_dir)
            error('%s resolves inside EEGLAB: %s', other_fns{k}, h);
        end
    end
end
