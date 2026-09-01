% s00_inventory.m  --  RUN THIS FIRST. (v2: fast, single-subject, no full
% pop_importbids call -- ds003702's recordings are long enough that importing
% all 50 subjects up front (as Pilot 1's s00 did for ds004117) can take a very
% long time with no console output in between. This version gets you the same
% answers in seconds by (a) reading channels.tsv / events.tsv directly with
% readtable -- exactly the format-agnostic approach already used in
% s00b_check_scaling.m's Method 1 -- and (b) loading only ONE subject's raw
% data via pop_loadbv for srate/channel/PSD checks.
%
% Nothing downstream (s01 onward) is valid until:
%   - cfg.ev_lock and cfg.maint_window are set from real trial timing
%   - cfg.cond_field / cfg.valid_field and their label values are correct
%   - cfg.roi_frontal / cfg.roi_parietal are confirmed against real channel labels
%   - cfg.resample_hz and cfg.line_hz are confirmed, not assumed

clear; clc;
cfg = config();

check_sub = 'sub-01';   % change to inspect a different subject

%% ---- PART A: direct TSV reads, no EEGLAB needed, seconds not minutes ----
sub_eeg_dir = fullfile(cfg.bids_dir, check_sub, 'eeg');
fprintf('\n========== ds003702 INVENTORY (%s) ==========\n', check_sub);
fprintf('Looking in: %s\n', sub_eeg_dir);

if ~exist(sub_eeg_dir, 'dir')
    error('Folder not found: %s -- check cfg.bids_dir in config.m.', sub_eeg_dir);
end

% --- channels.tsv ---
chanFile = dir(fullfile(sub_eeg_dir, '*_channels.tsv'));
if isempty(chanFile)
    warning('No *_channels.tsv found in %s', sub_eeg_dir);
else
    Tchan = readtable(fullfile(chanFile(1).folder, chanFile(1).name), ...
        'FileType','text','Delimiter','\t');
    fprintf('\n--- channels.tsv (%s) ---\n', chanFile(1).name);
    fprintf('Columns: %s\n', strjoin(Tchan.Properties.VariableNames, ', '));
    fprintf('n_channels = %d\n', height(Tchan));
    disp(Tchan.name');   % channel label list
    if ismember('units', Tchan.Properties.VariableNames)
        fprintf('Units found: %s\n', strjoin(unique(string(Tchan.units)), ', '));
    end
    if ismember('type', Tchan.Properties.VariableNames)
        fprintf('Channel types: %s\n', strjoin(unique(string(Tchan.type)), ', '));
    end
end

% --- events.tsv (THE key file for cfg.ev_lock / cond_field / valid_field / role_field) ---
evFile = dir(fullfile(sub_eeg_dir, '*_events.tsv'));
if isempty(evFile)
    warning('No *_events.tsv found in %s -- event info will only come from .vmrk (raw markers).', sub_eeg_dir);
else
    Tev = readtable(fullfile(evFile(1).folder, evFile(1).name), ...
        'FileType','text','Delimiter','\t');
    fprintf('\n--- events.tsv (%s) ---\n', evFile(1).name);
    fprintf('Columns: %s\n', strjoin(Tev.Properties.VariableNames, ', '));
    fprintf('n_events = %d\n', height(Tev));
    fprintf('\nFirst 20 rows:\n');
    disp(Tev(1:min(20,height(Tev)), :));

    fprintf('\nUnique values per column (up to 15 each):\n');
    for c = 1:numel(Tev.Properties.VariableNames)
        colname = Tev.Properties.VariableNames{c};
        col = Tev.(colname);
        try
            u = unique(col);
            if iscell(u) || isstring(u)
                u = u(1:min(15,numel(u)));
                fprintf('  %-20s : %s\n', colname, strjoin(string(u), ', '));
            elseif isnumeric(u)
                fprintf('  %-20s : numeric, range [%.4g, %.4g], %d unique values\n', ...
                    colname, min(u), max(u), numel(u));
            end
        catch
            fprintf('  %-20s : (could not summarize)\n', colname);
        end
    end

    % --- latency-gap hint using onset column (seconds, BIDS standard) ---
    if ismember('onset', Tev.Properties.VariableNames)
        fprintf('\n--- Time gaps between consecutive events (first 20 rows, ms) ---\n');
        onsets = Tev.onset(1:min(20,height(Tev)));
        for i = 2:numel(onsets)
            fprintf('  row %2d -> %2d : %8.1f ms\n', i-1, i, (onsets(i)-onsets(i-1))*1000);
        end
    end
end

%% ---- PART B: single-subject raw load for srate / channels / line-noise PSD ----
vhdrFile = dir(fullfile(sub_eeg_dir, '*_eeg.vhdr'));
if isempty(vhdrFile)
    warning('No *_eeg.vhdr found in %s -- skipping raw load.', sub_eeg_dir);
else
    fprintf('\n--- Loading raw BrainVision header/data for %s (single subject only) ---\n', check_sub);
    addpath(cfg.eeglab_dir); eeglab nogui;
    EEG = pop_loadbv(vhdrFile(1).folder, vhdrFile(1).name);

    fprintf('Sampling rate : %g Hz\n', EEG.srate);
    fprintf('n_channels    : %d\n', EEG.nbchan);
    fprintf('Duration      : %.1f s\n', EEG.pnts/EEG.srate);

    fprintf('\n--- RAW MARKER TYPES (.vmrk, generic BrainVision codes -- NOT the semantic\n');
    fprintf('    BIDS labels; use events.tsv above for real condition/validity/accuracy names) ---\n');
    if isfield(EEG,'event') && ~isempty(EEG.event)
        mtypes = unique({EEG.event.type});
        for i = 1:numel(mtypes)
            fprintf('  %-24s n=%d\n', mtypes{i}, sum(strcmp({EEG.event.type}, mtypes{i})));
        end
    end

    % raw amplitude sanity check (same spirit as s00b's scaling check)
    ch1 = double(EEG.data(1,:));
    fprintf('\n--- Raw amplitude check, channel 1 (%s) ---\n', EEG.chanlocs(1).labels);
    fprintf('  std = %.4g | range = %.4g\n', std(ch1), max(ch1)-min(ch1));

    % line-noise PSD check
    fprintf('\n--- Line-noise PSD check, channel 1 ---\n');
    try
        [pw, fr] = pwelch(ch1, [], [], [], EEG.srate);
        p50 = interp1(fr, pw, 50, 'nearest');
        p60 = interp1(fr, pw, 60, 'nearest');
        fprintf('  PSD power near 50 Hz = %.4g\n', p50);
        fprintf('  PSD power near 60 Hz = %.4g\n', p60);
        fprintf('  -> set cfg.line_hz to whichever shows the clear peak.\n');
    catch ME
        fprintf('  [WARN] pwelch check failed: %s\n', ME.message);
    end
end

fprintf(['\nACTION: set in config.m ->\n' ...
'  cfg.ev_lock                 = <event label to time-lock the connectivity window to>\n' ...
'  cfg.maint_window             = <[start end] s rel. to cfg.ev_lock, from the onset-gap>\n' ...
'                                   <printout above -- must land inside the delay period,>\n' ...
'                                   <before the probe>\n' ...
'  cfg.cond_field / cond_social / cond_nonsocial   = <social vs non-social cue field+values>\n' ...
'  cfg.valid_field / valid_valid / valid_invalid   = <cue-validity field+values, if used>\n' ...
'  cfg.role_field / role_feedback_correct / _incorrect / role_bad = <accuracy/trial-status field>\n' ...
'  cfg.resample_hz              = <pick a target rate <= the native rate seen above>\n' ...
'  cfg.line_hz                  = <50 or 60, from the PSD check above>\n' ...
'Also verify roi_frontal / roi_parietal labels exist in the channels.tsv list above.\n' ...
'\nIMPORTANT: repeat this same inventory on at least one more subject (change\n' ...
'check_sub above) before trusting any of the above as uniform across the dataset --\n' ...
'ds004117''s sub-022/sub-023 surprise came from checking only one subject first.\n' ...
'\nNOTE: this version deliberately avoids pop_importbids on the full dataset for\n' ...
'the inventory step, since importing all 50 subjects'' continuous data up front can\n' ...
'run for a very long time with no console output. s01_preprocess.m will still need\n' ...
'a decision on whether to keep using pop_importbids (with cfg.subset_n for a fast\n' ...
'first pass) or switch to a per-subject pop_loadbv + events.tsv-merge loop -- worth\n' ...
'deciding once we see how slow the full pop_importbids actually is on this dataset.\n']);
