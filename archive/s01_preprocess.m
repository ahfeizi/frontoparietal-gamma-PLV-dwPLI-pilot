% s01_preprocess.m  --  Per-SUBJECT preprocessing (Pilot 2 / ds003702).
% Unlike ds004117 (many run-datasets per subject via pop_importbids), each
% subject here has exactly ONE raw BrainVision file, so we load directly
% with pop_loadbv per subject -- much faster than importing the whole
% BIDS tree, and avoids the multi-minute-per-subject stall we hit earlier.
%
% Steps: load -> drop M1/M2/EOG/BIP1-24 -> filter [hp_hz lp_hz] -> line
% noise cleaning -> clean_rawdata (ASR + bad-channel detection) -> ICA ->
% ICLabel, with the gamma-specific Muscle-IC dominant-channel check ->
% interpolate flagged channels -> average reference -> save continuous
% clean .set + manifest (kept 1:1 key<->subject here, but manifest format
% matches Pilot 1's so downstream scripts don't need to special-case it).

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
if ~exist(cfg.deriv_dir,'dir'), mkdir(cfg.deriv_dir); end

sub_dirs = dir(fullfile(cfg.bids_dir, 'sub-*'));
sub_dirs = sub_dirs([sub_dirs.isdir]);
nAll = numel(sub_dirs);
nDo  = nAll; if ~isempty(cfg.subset_n), nDo = min(cfg.subset_n, nAll); end
fprintf('Subjects available: %d | processing: %d\n', nAll, nDo);

manifest_file = fullfile(cfg.deriv_dir,'manifest.mat');
failures_file = fullfile(cfg.deriv_dir,'preprocess_failures.mat');

% ---- RESUME SUPPORT: load prior progress if this run was interrupted ----
if exist(manifest_file, 'file')
    L = load(manifest_file); manifest = L.manifest;
    fprintf('Resuming: found existing manifest with %d completed subjects.\n', height(manifest));
else
    manifest = table('Size',[0 3], 'VariableTypes',{'string','string','string'}, ...
                     'VariableNames',{'key','subject','clean_file'});
end
if exist(failures_file, 'file')
    L2 = load(failures_file); failures = L2.failures;
else
    failures = table('Size',[0 3], 'VariableTypes',{'string','string','string'}, ...
                     'VariableNames',{'subject','step','error_message'});
end

for s = 1:nDo
    subj = sub_dirs(s).name;   % e.g. 'sub-01'
    key  = subj;               % 1:1 here (one file per subject)

    % ---- skip if already done (either a successful clean.set exists, or
    %      this subject already logged as a known failure) ----
    clean_file = fullfile(cfg.deriv_dir, [key '_clean.set']);
    already_done = any(manifest.subject == string(subj));
    already_failed = any(failures.subject == string(subj));
    if already_done || already_failed
        fprintf('=== %s  (%d/%d) -- already processed, skipping ===\n', key, s, nDo);
        continue;
    end
    if exist(clean_file, 'file') && ~already_done
        % .set exists on disk but manifest doesn't know about it (e.g. crash
        % happened right after saving but before this run's manifest was
        % updated) -- register it rather than reprocessing from scratch.
        fprintf('=== %s  (%d/%d) -- clean.set found on disk, registering without reprocessing ===\n', key, s, nDo);
        manifest = [manifest; {string(key), string(subj), string([key '_clean.set'])}]; %#ok<AGROW>
        save(manifest_file, 'manifest');
        writetable(manifest, fullfile(cfg.deriv_dir,'manifest.csv'));
        continue;
    end

    fprintf('\n=== %s  (%d/%d) ===\n', key, s, nDo);

    eeg_dir = fullfile(cfg.bids_dir, subj, 'eeg');
    vhdrFile = dir(fullfile(eeg_dir, '*_eeg.vhdr'));
    if isempty(vhdrFile)
        warning('%s: no .vhdr found in %s, skipping.', subj, eeg_dir);
        failures = [failures; {string(subj), "no_vhdr_found", ""}]; %#ok<AGROW>
        save(failures_file, 'failures');
        writetable(failures, fullfile(cfg.deriv_dir,'preprocess_failures.csv'));
        continue;
    end

    try

    EEG = pop_loadbv(vhdrFile(1).folder, vhdrFile(1).name);
    EEG.subject = subj;

    % ---- drop non-scalp channels (CONFIRMED list, from authors' script) ----
    keepChans = setdiff({EEG.chanlocs.labels}, cfg.exclude_channels, 'stable');
    EEG = pop_select(EEG, 'channel', keepChans);
    fprintf('  kept %d channels after excluding M1/M2/EOG/BIP1-24\n', EEG.nbchan);

    % ---- assign channel locations (pop_loadbv only reads labels, not
    %      positions -- ICLabel's topoplotFast-based feature extraction
    %      needs real theta/radius/X/Y/Z, which is why the earlier run
    %      failed inside pop_iclabel). Labels follow standard 10-10/10-5
    %      naming, so a standard lookup via the dipfit plugin works. ----
    elcCandidates = dir(fullfile(cfg.eeglab_dir, 'plugins', 'dipfit*', '**', 'standard-10-5-cap385.elc'));
    if isempty(elcCandidates)
        elcCandidates = dir(fullfile(cfg.eeglab_dir, 'plugins', 'dipfit*', '**', 'standard_1005.elc'));
    end
    if isempty(elcCandidates)
        error(['Could not find a standard channel location file under ' ...
            '%s/plugins/dipfit*. Verify the dipfit plugin is installed, ' ...
            'or set a location file path manually.'], cfg.eeglab_dir);
    end
    locfile = fullfile(elcCandidates(1).folder, elcCandidates(1).name);
    EEG = pop_chanedit(EEG, 'lookup', locfile);
    nMissingLoc = sum(cellfun(@isempty, {EEG.chanlocs.X}));
    if nMissingLoc > 0
        warning('%s: %d channel(s) got no location from the lookup -- check label spelling against %s.', ...
            subj, nMissingLoc, locfile);
    end

    % ---- resample (no-op if already at cfg.resample_hz -- see config.m note
    %      about trialinfo sample-index alignment) ----
    if EEG.srate ~= cfg.resample_hz
        warning(['%s: native srate %g Hz differs from cfg.resample_hz %g Hz. ' ...
            'Resampling will invalidate trialinfo sample indices unless you ' ...
            'scale them by (cfg.resample_hz/%g) in s02_epoch.m. Proceeding ' ...
            'with resample as configured -- make sure s02 accounts for this.'], ...
            subj, EEG.srate, cfg.resample_hz, EEG.srate);
        EEG = pop_resample(EEG, cfg.resample_hz);
    end

    % ---- filter (widened for gamma) ----
    EEG = pop_eegfiltnew(EEG, 'locutoff', cfg.hp_hz);
    EEG = pop_eegfiltnew(EEG, 'hicutoff', cfg.lp_hz);
    if ~isempty(cfg.line_hz)
        % NOTE: pop_cleanline requires the 'cleanline' EEGLAB plugin, which
        % is not installed on this machine (confirmed absent from the
        % eeglab nogui plugin list). Using a notch (band-stop) filter via
        % pop_eegfiltnew instead, which only needs the already-installed
        % firfilt plugin. Slightly less precise than adaptive cleanline,
        % but sufficient here -- install the cleanline plugin and swap
        % this back if you want the adaptive version.
        EEG = pop_eegfiltnew(EEG, 'locutoff', cfg.line_hz-1, ...
            'hicutoff', cfg.line_hz+1, 'revfilt', 1);
    end

    % ---- bad-channel / burst detection (same params as Pilot 1) ----
    EEG.urchanlocs = EEG.chanlocs;
    EEG = pop_clean_rawdata(EEG, 'FlatlineCriterion',5, 'ChannelCriterion',0.8, ...
        'LineNoiseCriterion',4, 'Highpass','off', 'BurstCriterion',20, ...
        'WindowCriterion','off', 'BurstRejection','off', 'Distance','Euclidian');

    % ---- ICA ----
    dataRank = sum(eig(cov(double(EEG.data'))) > 1e-7);
    % Switched from runica (extended infomax) to picard: several-fold faster
    % on this hardware (i7 4th gen, 12GB RAM) with equivalent decomposition
    % quality, and applied uniformly to ALL subjects (including the first 22
    % reprocessed after the runica run) so the ICA algorithm is consistent
    % across the full dataset -- not mixed between subjects.
    EEG = pop_runica(EEG, 'icatype','picard', 'pca', dataRank, 'mode','standard');
    EEG = pop_iclabel(EEG, 'default');

    % ICLabel column order: Brain Muscle Eye Heart LineNoise ChannelNoise Other
    probs = EEG.etc.ic_classification.ICLabel.classifications;
    muscle_prob = probs(:,2);
    eye_prob    = probs(:,3);
    heart_prob  = probs(:,4);
    line_prob   = probs(:,5);
    chan_prob   = probs(:,6);

    pth = cfg.iclabel.reject_prob;
    reject_mask = (eye_prob>=pth(1) & eye_prob<=pth(2)) ...
                | (heart_prob>=pth(1) & heart_prob<=pth(2)) ...
                | (line_prob>=pth(1) & line_prob<=pth(2)) ...
                | (chan_prob>=pth(1) & chan_prob<=pth(2));
    remove_ics = find(reject_mask)';   % non-Muscle artifact classes: remove outright

    % ---- gamma-specific Muscle-IC handling: dominant-channel check ----
    muscle_ics = find(muscle_prob >= cfg.muscle_ic.prob_thresh);
    interp_channels = {};
    for ic = muscle_ics(:)'
        topo = abs(EEG.icawinv(:, ic));
        z = (topo - mean(topo)) / std(topo);
        [sorted_topo, sort_idx] = sort(topo, 'descend');
        ratio = sorted_topo(1) / max(sorted_topo(2), eps);
        [maxz, maxz_idx] = max(z);

        if maxz >= cfg.muscle_ic.zthresh && ratio >= cfg.muscle_ic.ratio_thresh
            bad_chan = EEG.chanlocs(sort_idx(1)).labels;
            fprintf('  IC%d: Muscle-flagged but dominant-channel pattern (z=%.2f, ratio=%.2f) -> channel %s marked for interpolation, component KEPT\n', ...
                ic, maxz, ratio, bad_chan);
            interp_channels{end+1} = bad_chan; %#ok<SAGROW>
        else
            fprintf('  IC%d: Muscle-flagged, diffuse multichannel pattern (z=%.2f, ratio=%.2f) -> true EMG, component REMOVED\n', ...
                ic, maxz, ratio);
            remove_ics(end+1) = ic; %#ok<AGROW>
        end
    end
    remove_ics = unique(remove_ics);

    fprintf('  Removing %d components total (%d Muscle-true, %d other artifact classes)\n', ...
        numel(remove_ics), sum(ismember(remove_ics, muscle_ics)), ...
        numel(remove_ics) - sum(ismember(remove_ics, muscle_ics)));
    EEG = pop_subcomp(EEG, remove_ics, 0);

    % ---- interpolate: channels dropped by clean_rawdata + dominant-channel
    %      Muscle-IC flags, all against the ORIGINAL (pre-clean_rawdata)
    %      channel locations ----
    interp_channels = unique(interp_channels);
    if ~isempty(interp_channels)
        fprintf('  Interpolating %d dominant-channel Muscle-IC channels: %s\n', ...
            numel(interp_channels), strjoin(interp_channels, ', '));
    end
    EEG = pop_interp(EEG, EEG.urchanlocs, 'spherical');

    % ---- average reference ----
    EEG = pop_reref(EEG, []);

    pop_saveset(EEG, 'filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    manifest = [manifest; {string(key), string(subj), string([key '_clean.set'])}]; %#ok<AGROW>
    save(manifest_file, 'manifest');
    writetable(manifest, fullfile(cfg.deriv_dir,'manifest.csv'));

    clear EEG;

    catch ME
        fprintf('  [FAILED] %s: %s\n', subj, ME.message);
        failures = [failures; {string(subj), string(ME.identifier), string(ME.message)}]; %#ok<AGROW>
        save(failures_file, 'failures');
        writetable(failures, fullfile(cfg.deriv_dir,'preprocess_failures.csv'));
        clear EEG;
        continue;
    end
end

writetable(manifest, fullfile(cfg.deriv_dir,'manifest.csv'));
save(fullfile(cfg.deriv_dir,'manifest.mat'),'manifest');
if height(failures) > 0
    writetable(failures, fullfile(cfg.deriv_dir,'preprocess_failures.csv'));
    fprintf('\n%d subject(s) FAILED and were skipped -- see preprocess_failures.csv:\n', height(failures));
    disp(failures);
end
fprintf('\nPreprocessing complete: %d subjects succeeded, %d failed.\n', height(manifest), height(failures));
