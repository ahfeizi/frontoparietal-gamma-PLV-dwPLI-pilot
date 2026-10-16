% s23_imc_diagnostic.m -- EXPLORATORY. Per-subject estimator diagnostic for
% ImC and dwPLI: how much each depends on averaging the cross-spectrum
% jointly over samples x trials versus across trials at each sample.
% gamma_high only; decompositions 'hilbert' and 'wavelet' (same bins as the
% battery); references 'average' and 'csd'; all 4 cue x validity cells;
% all 64 frontal x parietal pairs. Values are signed (no abs).
%
% With X = zi.*conj(zj) as a [samples x trials] matrix, per pair:
%   imc_joint   : imag(Sxy/sqrt(Sxx*Syy)), spectra averaged jointly over
%                 samples and trials -- exactly util_conn_pairs.m (ImC)
%   imc_trial   : at each sample imag(mean_k X)/sqrt(mean_k|zi|^2 * mean_k|zj|^2)
%                 (k = trials), then mean over samples
%   dwpli_trial : dwPLI across trials at each sample, mean over samples --
%                 exactly util_conn_pairs.m (dwPLI)
%   dwpli_joint : the same debiased estimator, with all samples x trials
%                 treated as one set of observations
% Per ROI channel (16): power = mean |z|^2 over samples and trials.
% Wavelet: each measure (and power) is averaged over the bins with >= 5
% trials, sum / n_used, as s21_battery.m.
%
% Data loading, trial selection, window, CSD, Hilbert and wavelet code are
% copied from s21_battery.m (local functions cannot be shared between
% scripts); imc_joint and dwpli_trial come from util_conn_pairs.m.
% s24_imc_diagnostic_check.m verifies them against battery_persubject.csv.
%
% Outputs (cfg.deriv_dir, appended one subject at a time):
%   imc_diag_persubject.csv  subject,condition,decomposition,reference,pair,measure,value
%   imc_diag_power.csv       subject,condition,decomposition,reference,channel,power
%   imc_diag_done_subjects.txt  completed subjects (resume bookkeeping)
% RESUME: as s21 -- completed subjects are skipped; rows of a subject that
% was interrupted mid-write are removed on restart. To start from scratch,
% delete the two CSVs and imc_diag_done_subjects.txt.

clear; clc;
repo = fileparts(fileparts(mfilename('fullpath')));   % repository root (this script is in battery/)
addpath(repo, fullfile(repo, 'battery'));             % config.m, util_bandpass.m, util_chanidx.m; util_conn_pairs.m
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

srate = cfg.resample_hz;
bn    = 'gamma_high';
band  = cfg.battery.bands.(bn);
freqs = cfg.battery.wavelet_freqs.(bn);
n_cycles = cfg.battery.n_cycles;
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
refs = {'average','csd'};
measNames = {'imc_joint','imc_trial','dwpli_trial','dwpli_joint'};

f_main  = fullfile(cfg.deriv_dir,'imc_diag_persubject.csv');
f_power = fullfile(cfg.deriv_dir,'imc_diag_power.csv');
f_done  = fullfile(cfg.deriv_dir,'imc_diag_done_subjects.txt');
done = local_resume(f_done, {f_main, f_power}, { ...
    'subject,condition,decomposition,reference,pair,measure,value', ...
    'subject,condition,decomposition,reference,channel,power'});

for su = 1:height(manifest)
    subj = char(manifest.subject(su));
    key  = char(manifest.key(su));
    if ismember(subj, done), fprintf('%s already complete, skipping\n', subj); continue; end
    t0 = tic;

    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    T = load(fullfile(cfg.deriv_dir, [key '_trials.mat'])); trials = T.trials;
    good_trials = trials([trials.good]);
    if isempty(good_trials)
        warning('%s: no good trials, skipping subject entirely.', subj);
        continue;
    end

    have = {EEG.chanlocs.labels};
    labF = local_present(have, cfg.roi_frontal);
    labP = local_present(have, cfg.roi_parietal);
    if numel(labF) < numel(cfg.roi_frontal) || numel(labP) < numel(cfg.roi_parietal)
        warning('%s: only %d/%d frontal and %d/%d parietal ROI channels found -- check labels.', ...
            subj, numel(labF), numel(cfg.roi_frontal), numel(labP), numel(cfg.roi_parietal));
    end
    roiLabs = [labF labP];

    [csd_data, csd_labels] = local_csd(EEG);
    iAvg = util_chanidx(EEG, roiLabs);
    iCsd = util_chanidx(local_labels2eeg(csd_labels), roiLabs);
    if numel(iCsd) ~= numel(roiLabs)
        error('%s: ROI channel(s) missing in the CSD output.', subj);
    end

    runs = struct('srate', srate, 'sel', {cell(1,nC)}, ...
        'average', double(EEG.data(iAvg,:)), 'csd', double(csd_data(iCsd,:)));
    for c = 1:nC
        runs.sel{c} = good_trials([good_trials.cue] == condCue(c) & [good_trials.valid] == condVal(c));
    end
    clear EEG csd_data

    nF = numel(labF); nP = numel(labP); np = nF*nP; nR = numel(roiLabs);
    fr = 1:nF; pa = nF + (1:nP);            % indices into the ROI-row segments
    pairLabs = cell(1,np);
    for i = 1:nF, for j = 1:nP, pairLabs{(i-1)*nP+j} = sprintf('%s-%s', labF{i}, labP{j}); end, end

    buf_main = {}; buf_power = {};
    for r = 1:numel(refs)
        rf = refs{r};

        % ---------------- band-pass + Hilbert (as s21 / s03) ----------------
        seg = local_pool(runs, rf, cfg.maint_window, ...
            @(X, sr) local_hilbert(X, sr, band, cfg.filt_order_factor));
        for c = 1:nC
            kept = local_ntr(seg{c});
            if kept < 5
                fprintf('  %s %s %s hilbert: only %d in-bounds segments, skipping cell\n', subj, rf, conds{c}, kept);
                continue;
            end
            [Mv, pw] = local_measures(seg{c}, fr, pa);
            buf_main{end+1}  = local_rows_main(subj, conds{c}, 'hilbert', rf, pairLabs, measNames, Mv); %#ok<SAGROW>
            buf_power{end+1} = local_rows_power(subj, conds{c}, 'hilbert', rf, roiLabs, pw); %#ok<SAGROW>
        end

        % ---------------- Morlet wavelet (as s21) ----------------
        nf = numel(freqs);
        W = cell(1,nC);
        for c = 1:nC
            W{c} = struct('ok', false(1,nf), 'M', nan(nf, np, numel(measNames)), 'pw', nan(nf, nR));
        end
        for f = 1:nf
            fc = freqs(f);
            seg = local_pool(runs, rf, cfg.maint_window, ...
                @(X, sr) local_morlet(X, sr, fc, n_cycles));
            for c = 1:nC
                if local_ntr(seg{c}) < 5, continue; end
                [Mv, pw] = local_measures(seg{c}, fr, pa);
                W{c}.M(f,:,:) = reshape(Mv, 1, np, numel(measNames));
                W{c}.pw(f,:)  = pw;
                W{c}.ok(f) = true;
            end
        end
        for c = 1:nC
            used = find(W{c}.ok);               % bins with >= 5 trials, as s21
            nu = numel(used);
            if nu == 0
                fprintf('  %s %s %s wavelet: no bin with >= 5 trials, skipping cell\n', subj, rf, conds{c});
                continue;
            end
            % per-pair (per-channel) average over bins, sum / n_used as s21 (no omitnan)
            Mv = reshape(sum(W{c}.M(used,:,:),1) / nu, np, numel(measNames));
            pw = sum(W{c}.pw(used,:),1) / nu;
            buf_main{end+1}  = local_rows_main(subj, conds{c}, 'wavelet', rf, pairLabs, measNames, Mv); %#ok<SAGROW>
            buf_power{end+1} = local_rows_power(subj, conds{c}, 'wavelet', rf, roiLabs, pw); %#ok<SAGROW>
        end
    end

    % ---- append this subject, then mark it complete ----
    local_append(f_main,  buf_main);
    local_append(f_power, buf_power);
    local_append(f_done,  {sprintf('%s\n', subj)});
    fprintf('imc diagnostic done %s (%.1f min)\n', subj, toc(t0)/60);
end
fprintf('\nIMC diagnostic complete. Outputs in %s\n', cfg.deriv_dir);

% ======================= local functions =======================
function [Mv, pw] = local_measures(seg, fr, pa)
    % Mv: [np x 4] = imc_joint, imc_trial, dwpli_trial, dwpli_joint
    % pw: [1 x nROI] mean |z|^2 over samples and trials, per ROI row of seg
    [~, dwpli_trial, imc_joint] = util_conn_pairs(seg, fr, pa);
    wlen = size(seg,2); kept = size(seg,3);
    np = numel(fr) * numel(pa);
    imc_trial = zeros(1,np); dwpli_joint = zeros(1,np); p = 0;
    for i = 1:numel(fr)
        zi = reshape(seg(fr(i),:,:), wlen, kept);
        for j = 1:numel(pa)
            p = p + 1;
            zj = reshape(seg(pa(j),:,:), wlen, kept);
            X = zi .* conj(zj);                               % [samples x trials]
            % ImC per sample across trials, then mean over samples
            Sxy_t = mean(X, 2);
            Sxx_t = mean(abs(zi).^2, 2);
            Syy_t = mean(abs(zj).^2, 2);
            imc_trial(p) = mean(imag(Sxy_t) ./ sqrt(Sxx_t .* Syy_t));
            % dwPLI with all samples x trials as one set of observations
            im = imag(X(:));
            s_im = sum(im); s_abs = sum(abs(im)); s_sq = sum(im.^2);
            dwpli_joint(p) = (s_im^2 - s_sq) / (s_abs^2 - s_sq);
        end
    end
    Mv = [imc_joint(:), imc_trial(:), dwpli_trial(:), dwpli_joint(:)];
    pw = mean(reshape(abs(seg).^2, size(seg,1), []), 2).';
end

function [csd_data, csd_labels] = local_csd(EEG)
    % ---- CSD, as s18_csd_sensitivity.m / s21_battery.m ----
    ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
    cfg_csd = []; cfg_csd.method = 'spline';
    ft_csd = ft_scalpcurrentdensity(cfg_csd, ft_raw);
    csd_data = ft_csd.trial{1};             % [nchan x ntime], CSD-transformed continuous data
    csd_labels = ft_csd.label;
end

function Z = local_hilbert(X, srate, band, order_factor)
    % as s03_connectivity.m
    Xf = util_bandpass(X, srate, band, order_factor);
    Z  = hilbert(Xf.').';
end

function Z = local_morlet(data, srate, fc, n_cycles)
    % complex Morlet convolution, as s21_battery.m (from frontotemporal s16)
    sigma_t = n_cycles / (2*pi*fc);
    t = -3*sigma_t : 1/srate : 3*sigma_t;
    wavelet = exp(2i*pi*fc.*t) .* exp(-t.^2/(2*sigma_t^2));
    wavelet = wavelet / sqrt(sum(abs(wavelet).^2));   % unit energy
    Z = zeros(size(data));
    for ch = 1:size(data,1)
        Z(ch,:) = conv(data(ch,:), wavelet, 'same');
    end
end

function seg = local_pool(runs, rf, maint_window, transform)
    % Maintenance-window segments of the selected trials, per condition
    % (as s03: trials out of the recording bounds are dropped).
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

function E = local_labels2eeg(labels)
    % minimal EEG-like struct so util_chanidx can index FieldTrip labels
    E.chanlocs = struct('labels', reshape(labels, 1, []));
end

function s = local_rows_main(subj, cond, dec, rf, pairLabs, measNames, Mv)
    np = numel(pairLabs); s = '';
    for p = 1:np
        C = [repmat({subj, cond, dec, rf, pairLabs{p}}, numel(measNames), 1), ...
            measNames(:), num2cell(Mv(p,:)')]';
        s = [s sprintf('%s,%s,%s,%s,%s,%s,%.17g\n', C{:})]; %#ok<AGROW>
    end
end

function s = local_rows_power(subj, cond, dec, rf, roiLabs, pw)
    nR = numel(roiLabs);
    C = [repmat({subj, cond, dec, rf}, nR, 1), roiLabs(:), num2cell(pw(:))]';
    s = sprintf('%s,%s,%s,%s,%s,%.17g\n', C{:});
end

function local_append(f, buf)
    fid = fopen(f, 'a');
    if fid < 0, error('cannot open %s for appending', f); end
    fprintf(fid, '%s', buf{:});
    fclose(fid);
end

function done = local_resume(f_done, files, headers)
    done = {};
    if isfile(f_done)
        done = strtrim(splitlines(fileread(f_done)));
        done = reshape(done(~cellfun(@isempty, done)), 1, []);
    end
    for i = 1:numel(files)
        if ~isfile(files{i})
            if ~isempty(done)
                error(['%s is missing but %s lists completed subjects. Delete %s ' ...
                    '(and the other imc_diag CSV) to restart from scratch.'], files{i}, f_done, f_done);
            end
            fid = fopen(files{i}, 'w'); fprintf(fid, '%s\n', headers{i}); fclose(fid);
            continue;
        end
        L = splitlines(fileread(files{i}));
        L = L(~cellfun(@isempty, L));
        body = L(2:end);
        subj = regexp(body, '^[^,]*', 'match', 'once');
        keep = ismember(subj, done);
        if all(keep), continue; end
        fprintf('Resume: removing %d rows of unfinished subject(s) from %s\n', sum(~keep), files{i});
        fid = fopen(files{i}, 'w'); fprintf(fid, '%s\n', L{1}, body{keep}); fclose(fid);
    end
    if ~isempty(done), fprintf('Resume: %d subject(s) already complete.\n', numel(done)); end
end
