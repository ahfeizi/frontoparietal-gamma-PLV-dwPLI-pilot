% s21_battery.m -- Convergent-validation battery, per-subject values only
% (NO statistics). Every combination of
%   measure       : PLV, dwPLI, ImC (signed imag(Sxy/sqrt(Sxx*Syy)), as s13)
%   decomposition : 'hilbert' (band-pass + Hilbert, as s03_connectivity.m)
%                   'wavelet' (complex Morlet, 5 cycles, ported from the
%                              frontotemporal pilot's s16_wavelet_beta.m;
%                              per-bin connectivity averaged within band)
%   reference     : 'average' (existing _clean.set) or 'csd' (FieldTrip
%                   ft_scalpcurrentdensity 'spline', as s18_csd_sensitivity.m)
% for every band in cfg.battery.bands (gamma_low, gamma_high), every
% subject, the four cue x validity cells (avatar_valid, avatar_invalid,
% stick_valid, stick_invalid) and all 64 frontal x parietal pairs.
%
% Data loading, trial selection (good trials from <key>_trials.mat:
% accuracy + RT filtering per s02_epoch.m), maintenance window (locked to
% trialinfo column Maint) and the >= 5 in-bounds-trials rule follow
% s03_connectivity.m. One continuous recording per subject, so no run
% pooling. PLV/dwPLI/ImC come from util_conn_pairs.m, so Hilbert and
% wavelet paths use one estimator.
%
% Outputs (in cfg.deriv_dir, long format, appended one subject at a time):
%   battery_persubject.csv     subject,condition,band,measure,decomposition,
%                              reference,pair,value
%   battery_phasediff.csv      subject,condition,band,decomposition,reference,
%                              pair,mean_angle_deg,resultant_len
%                              (resultant = mean(exp(1i*dphi)) over samples x
%                              trials, as s13; for the wavelet path the per-bin
%                              complex vectors are averaged across the band's bins)
%   battery_phasediff_bins.csv subject,condition,band,reference,pair,freq_hz,
%                              mean_angle_deg,resultant_len (wavelet only, one
%                              row per centre-frequency bin)
%   battery_done_subjects.txt  completed subjects (resume bookkeeping)
%
% RESUME: a subject is written to all three CSVs and then to the done list.
% On restart, completed subjects are skipped and any rows of a subject that
% was interrupted mid-write are removed first. To start from scratch, delete
% the three CSVs and battery_done_subjects.txt.

clear; clc;
repo = fileparts(fileparts(mfilename('fullpath')));   % repository root (this script is in battery/)
addpath(repo, fullfile(repo, 'battery'));             % config.m, util_bandpass.m, util_chanidx.m; util_conn_pairs.m
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

srate = cfg.resample_hz;
bands = fieldnames(cfg.battery.bands);
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
refs  = {'average','csd'};
n_cycles = cfg.battery.n_cycles;
all_freqs = [];
for b = 1:numel(bands), all_freqs = [all_freqs cfg.battery.wavelet_freqs.(bands{b})]; end %#ok<AGROW>
all_freqs = unique(all_freqs);

f_main  = fullfile(cfg.deriv_dir,'battery_persubject.csv');
f_phase = fullfile(cfg.deriv_dir,'battery_phasediff.csv');
f_bins  = fullfile(cfg.deriv_dir,'battery_phasediff_bins.csv');
f_done  = fullfile(cfg.deriv_dir,'battery_done_subjects.txt');
done = local_resume(f_done, {f_main, f_phase, f_bins}, { ...
    'subject,condition,band,measure,decomposition,reference,pair,value', ...
    'subject,condition,band,decomposition,reference,pair,mean_angle_deg,resultant_len', ...
    'subject,condition,band,reference,pair,freq_hz,mean_angle_deg,resultant_len'});

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

    % one "run" per subject; same pooling code path as the frontotemporal battery
    runs = struct('srate', srate, 'sel', {cell(1,nC)}, ...
        'average', double(EEG.data(iAvg,:)), 'csd', double(csd_data(iCsd,:)));
    for c = 1:nC
        runs.sel{c} = good_trials([good_trials.cue] == condCue(c) & [good_trials.valid] == condVal(c));
    end
    clear EEG csd_data

    nF = numel(labF); nP = numel(labP); np = nF*nP;
    fr = 1:nF; pa = nF + (1:nP);            % indices into the ROI-row segments
    pairLabs = cell(1,np);
    for i = 1:nF, for j = 1:nP, pairLabs{(i-1)*nP+j} = sprintf('%s-%s', labF{i}, labP{j}); end, end

    buf_main = {}; buf_phase = {}; buf_bins = {};
    for r = 1:numel(refs)
        rf = refs{r};

        % ---------------- band-pass + Hilbert (as s03) ----------------
        for b = 1:numel(bands)
            bn = bands{b}; band = cfg.battery.bands.(bn);
            seg = local_pool(runs, rf, cfg.maint_window, ...
                @(X, sr) local_hilbert(X, sr, band, cfg.filt_order_factor));
            for c = 1:nC
                kept = local_ntr(seg{c});
                if kept < 5
                    fprintf('  %s %s %s %s hilbert: only %d in-bounds segments, skipping cell\n', ...
                        subj, rf, bn, conds{c}, kept);
                    continue;
                end
                [plv, dwpli, imc, mvec] = util_conn_pairs(seg{c}, fr, pa);
                buf_main{end+1}  = local_rows_main(subj, conds{c}, bn, 'hilbert', rf, pairLabs, plv, dwpli, imc); %#ok<SAGROW>
                buf_phase{end+1} = local_rows_phase(subj, conds{c}, bn, 'hilbert', rf, pairLabs, mvec); %#ok<SAGROW>
            end
        end

        % ---------------- Morlet wavelet (ported from frontotemporal s16) ----------------
        nf = numel(all_freqs);
        W = cell(1,nC);
        for c = 1:nC
            W{c} = struct('ok', false(1,nf), 'plv', nan(nf,np), 'dwpli', nan(nf,np), ...
                'imc', nan(nf,np), 'mvec', complex(nan(nf,np)));
        end
        for f = 1:nf
            fc = all_freqs(f);
            seg = local_pool(runs, rf, cfg.maint_window, ...
                @(X, sr) local_morlet(X, sr, fc, n_cycles));
            for c = 1:nC
                if local_ntr(seg{c}) < 5, continue; end
                [W{c}.plv(f,:), W{c}.dwpli(f,:), W{c}.imc(f,:), W{c}.mvec(f,:)] = ...
                    util_conn_pairs(seg{c}, fr, pa);
                W{c}.ok(f) = true;
            end
        end
        for b = 1:numel(bands)
            bn = bands{b};
            fb = find(ismember(all_freqs, cfg.battery.wavelet_freqs.(bn)));
            for c = 1:nC
                used = fb(W{c}.ok(fb));        % bins with >= 5 trials, as s16
                nu = numel(used);
                if nu == 0
                    fprintf('  %s %s %s %s wavelet: no bin with >= 5 trials, skipping cell\n', ...
                        subj, rf, bn, conds{c});
                    continue;
                end
                % per-pair average over bins, sum / n_used as s16 (no omitnan)
                plv   = sum(W{c}.plv(used,:),1)   / nu;
                dwpli = sum(W{c}.dwpli(used,:),1) / nu;
                imc   = sum(W{c}.imc(used,:),1)   / nu;
                mvec  = sum(W{c}.mvec(used,:),1)  / nu;
                buf_main{end+1}  = local_rows_main(subj, conds{c}, bn, 'wavelet', rf, pairLabs, plv, dwpli, imc); %#ok<SAGROW>
                buf_phase{end+1} = local_rows_phase(subj, conds{c}, bn, 'wavelet', rf, pairLabs, mvec); %#ok<SAGROW>
                for u = used
                    buf_bins{end+1} = local_rows_bins(subj, conds{c}, bn, rf, pairLabs, all_freqs(u), W{c}.mvec(u,:)); %#ok<SAGROW>
                end
            end
        end
    end

    % ---- append this subject, then mark it complete ----
    local_append(f_main,  buf_main);
    local_append(f_phase, buf_phase);
    local_append(f_bins,  buf_bins);
    local_append(f_done,  {sprintf('%s\n', subj)});
    fprintf('battery done %s (%.1f min)\n', subj, toc(t0)/60);
end
fprintf('\nBattery complete. Outputs in %s\n', cfg.deriv_dir);

% ======================= local functions =======================
function [csd_data, csd_labels] = local_csd(EEG)
    % ---- CSD, as s18_csd_sensitivity.m ----
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
    % complex Morlet convolution, lines as frontotemporal s16_wavelet_beta.m
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

function s = local_rows_main(subj, cond, band, dec, rf, pairLabs, plv, dwpli, imc)
    meas = {'PLV','dwPLI','ImC'}; vals = {plv, dwpli, imc};
    np = numel(pairLabs); s = '';
    for m = 1:numel(meas)
        C = [repmat({subj, cond, band, meas{m}, dec, rf}, np, 1), pairLabs(:), num2cell(vals{m}(:))]';
        s = [s sprintf('%s,%s,%s,%s,%s,%s,%s,%.17g\n', C{:})]; %#ok<AGROW>
    end
end

function s = local_rows_phase(subj, cond, band, dec, rf, pairLabs, mvec)
    np = numel(pairLabs);
    C = [repmat({subj, cond, band, dec, rf}, np, 1), pairLabs(:), ...
        num2cell(angle(mvec(:))*180/pi), num2cell(abs(mvec(:)))]';
    s = sprintf('%s,%s,%s,%s,%s,%s,%.17g,%.17g\n', C{:});
end

function s = local_rows_bins(subj, cond, band, rf, pairLabs, fc, mvec)
    np = numel(pairLabs);
    C = [repmat({subj, cond, band, rf}, np, 1), pairLabs(:), repmat({fc}, np, 1), ...
        num2cell(angle(mvec(:))*180/pi), num2cell(abs(mvec(:)))]';
    s = sprintf('%s,%s,%s,%s,%s,%.17g,%.17g,%.17g\n', C{:});
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
                    '(and the other battery CSVs) to restart from scratch.'], files{i}, f_done, f_done);
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
