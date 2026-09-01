% s03_connectivity.m  --  THE SCIENTIFIC CORE (subject-level), Pilot 2.
% For each subject, for each gamma sub-band (gamma_low, gamma_high) and each
% of the 2x2 cue (avatar/stick) x validity (valid/invalid) cells, extracts
% the maintenance-window analytic segments (cfg.maint_window, locked to
% trialinfo column Maint) from GOOD trials only (per s02_epoch.m: accuracy +
% RT filtering, boundary-crossing and out-of-bounds trials excluded), then
% computes:
%   dwPLI : debiased weighted phase-lag index (Vinck et al., 2011)  [PRIMARY]
%   PLV   : phase-locking value                                     [SECONDARY]
% over all frontal x parietal channel pairs, plus mean band power within
% each ROI separately (frontal_power, parietal_power) as a local-power
% control analogous to Pilot 1's frontal-midline power measure.
% One row per subject x band x cue x validity. N for statistics = subjects.
%
% Unlike Pilot 1 (which pooled trials across multiple per-subject run
% files), each subject here has exactly one continuous recording, so no
% cross-run pooling is needed.
%
% Requires util_bandpass.m and util_chanidx.m on the path (same utilities
% used by Pilot 1's s03_connectivity.m).

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;

M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

srate = cfg.resample_hz;
win = round(cfg.maint_window * srate);
wlen = win(2) - win(1) + 1;

bands   = fieldnames(cfg.bands);           % {'gamma_low','gamma_high'}
cueLabs = {'stick','avatar'};
cueVals = [cfg.ti.cue_stick, cfg.ti.cue_avatar];
valLabs = {'valid','invalid'};
valVals = [cfg.ti.val_valid, cfg.ti.val_invalid];

rows = {};

for su = 1:height(manifest)
    subj = char(manifest.subject(su));
    key  = char(manifest.key(su));
    fprintf('\n=== %s ===\n', subj);

    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    T = load(fullfile(cfg.deriv_dir, [key '_trials.mat'])); trials = T.trials;

    fr = util_chanidx(EEG, cfg.roi_frontal);
    pa = util_chanidx(EEG, cfg.roi_parietal);
    if numel(fr) < numel(cfg.roi_frontal)
        warning('%s: only %d/%d frontal ROI channels found -- check labels.', ...
            subj, numel(fr), numel(cfg.roi_frontal));
    end
    if numel(pa) < numel(cfg.roi_parietal)
        warning('%s: only %d/%d parietal ROI channels found -- check labels.', ...
            subj, numel(pa), numel(cfg.roi_parietal));
    end

    good_trials = trials([trials.good]);
    if isempty(good_trials)
        warning('%s: no good trials, skipping subject entirely.', subj);
        continue;
    end

    for b = 1:numel(bands)
        bandName = bands{b};
        Xf = util_bandpass(double(EEG.data), srate, cfg.bands.(bandName), cfg.filt_order_factor);
        Z  = hilbert(Xf.').';   % [nchan x ntime] analytic signal

        for c = 1:numel(cueLabs)
            for v = 1:numel(valLabs)
                sel = good_trials([good_trials.cue] == cueVals(c) & [good_trials.valid] == valVals(v));
                ntr = numel(sel);
                if ntr < 5
                    fprintf('  %s %s-%s: only %d good trials, skipping cell\n', ...
                        bandName, cueLabs{c}, valLabs{v}, ntr);
                    continue;
                end

                seg = nan(EEG.nbchan, wlen, ntr);
                kept = 0;
                for t = 1:ntr
                    a = sel(t).lock_sample + win(1);
                    z = a + wlen - 1;
                    if a < 1 || z > size(Z,2), continue; end
                    kept = kept + 1;
                    seg(:,:,kept) = Z(:, a:z);
                end
                seg = seg(:,:,1:kept);
                if kept < 5
                    fprintf('  %s %s-%s: only %d in-bounds segments, skipping cell\n', ...
                        bandName, cueLabs{c}, valLabs{v}, kept);
                    continue;
                end

                % ---- dwPLI / PLV over all frontal x parietal pairs ----
                np = numel(fr) * numel(pa);
                dwpli_pairs = zeros(1,np); plv_pairs = zeros(1,np); p = 0;
                for i = 1:numel(fr)
                    zi = reshape(seg(fr(i),:,:), wlen, kept);
                    for j = 1:numel(pa)
                        p = p + 1;
                        zj = reshape(seg(pa(j),:,:), wlen, kept);
                        X = zi .* conj(zj); im = imag(X);
                        s_im = sum(im,2); s_abs = sum(abs(im),2); s_sq = sum(im.^2,2);
                        dwpli_t = (s_im.^2 - s_sq) ./ (s_abs.^2 - s_sq);
                        plv_t   = abs(mean(X ./ abs(X), 2));
                        dwpli_pairs(p) = mean(dwpli_t, 'omitnan');
                        plv_pairs(p)   = mean(plv_t, 'omitnan');
                    end
                end

                % ---- local band power within each ROI (control measure) ----
                frontal_power  = mean(mean(abs(reshape(seg(fr,:,:), numel(fr), wlen*kept)).^2, 2));
                parietal_power = mean(mean(abs(reshape(seg(pa,:,:), numel(pa), wlen*kept)).^2, 2));

                rows(end+1,:) = {subj, bandName, cueLabs{c}, valLabs{v}, ...
                    mean(dwpli_pairs,'omitnan'), mean(plv_pairs,'omitnan'), ...
                    frontal_power, parietal_power, kept, {dwpli_pairs}, {plv_pairs}}; %#ok<SAGROW>

                fprintf('  %s %s-%s: n=%d trials | dwPLI=%.4f PLV=%.4f frontal_pow=%.3g parietal_pow=%.3g\n', ...
                    bandName, cueLabs{c}, valLabs{v}, kept, ...
                    mean(dwpli_pairs,'omitnan'), mean(plv_pairs,'omitnan'), frontal_power, parietal_power);
            end
        end
    end
end

R = cell2table(rows, 'VariableNames', {'subject','band','cue','valid', ...
    'dwpli','plv','frontal_power','parietal_power','ntrials','dwpli_pairs','plv_pairs'});
save(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R', 'cfg');
writetable(R(:,1:9), fullfile(cfg.deriv_dir,'conn_results.csv'));
fprintf('\nSaved conn_results: %d subjects.\n', numel(unique(R.subject)));
