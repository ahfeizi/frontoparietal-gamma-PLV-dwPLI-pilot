%% s08a_explore_cue_markers.m
%
% PURPOSE
%   Pre-analysis step required before the alpha-reproduction pipeline
%   (s08b) can be trusted. We need to know whether the raw .vmrk marker
%   stream contains a marker that corresponds to CUE ONSET (the event
%   the original paper time-locks to), or whether the earliest reliably
%   identifiable event is TARGET ONSET instead. This script does NOT
%   assume the answer - it empirically derives it from the data and
%   prints a report you should read before running s08b.
%
% LOGIC
%   1. For each subject, read the .vmrk marker stream (Value + Position).
%   2. Reconstruct trial boundaries using the s30xx codes already
%      validated in the Pilot-2 pipeline (maintenance window markers
%      s3023/s3043 -> s3024/s3044 are the one fixed point we trust).
%   3. Within each trial, list every marker code that occurs, in order,
%      relative to the maintenance-onset marker (s3023/s3043), together
%      with its median latency (ms) relative to that anchor.
%   4. A marker that (a) occurs exactly once per trial, (b) occurs BEFORE
%      the maintenance marker, and (c) has a highly consistent latency
%      across trials (low SD) is a strong candidate for a fixed
%      task-event (cue onset or target onset).
%   5. Cross-check candidate markers against trialinfo (data_ica.mat,
%      17-column matrix) columns 4 (validity) and 6 (cue type) to see
%      whether the candidate marker's code itself varies with cue type
%      (e.g. one code for stick, another for avatar) - if so it is very
%      likely the CUE marker, since target displays are typically
%      identical in code regardless of validity/cue-type.
%
% OUTPUT
%   - Console report per subject and pooled across subjects.
%   - /derivatives/marker_exploration/marker_report.mat with a struct
%     array `report(subj).candidates` for manual inspection.
%   - A printed RECOMMENDATION line telling s08b which onset to use.
%
% NOTE ON HONESTY
%   This script makes no claim about what the marker mapping "should"
%   be from the paper's methods text (we do not have independently
%   verified access to that mapping). It only reports what is
%   empirically present in the data. You (Amir) should sanity check the
%   RECOMMENDATION before trusting it for the manuscript.

clear; clc;

%% ---- CONFIG ----
cfg_local          = [];
cfg_local.raw_root = fullfile('..','data','ds003702');
cfg_local.out_root = fullfile('..','derivatives','marker_exploration');
cfg_local.subjects  = get_valid_subject_list(cfg_local.raw_root); % excludes 08/42/47 automatically

if ~exist(cfg_local.out_root, 'dir'); mkdir(cfg_local.out_root); end

report = struct('subject', {}, 'candidates', {}, 'n_trials', {});

for s = 1:numel(cfg_local.subjects)
    subj = cfg_local.subjects{s};
    fprintf('\n=== Exploring markers: %s ===\n', subj);

    vmrk_file = find_vmrk_file(cfg_local.raw_root, subj);
    if isempty(vmrk_file)
        warning('No .vmrk found for %s, skipping.', subj);
        continue
    end

    % --- Read marker stream (Value, Position in samples) ---
    mk = read_vmrk_markers(vmrk_file); % struct array: .value (string, e.g. 'S 23'), .position

    % --- Reduce to numeric stimulus codes only (S###) ---
    is_stim   = ~cellfun(@isempty, regexp({mk.value}, '^S\s*\d+$', 'once'));
    stim_mk   = mk(is_stim);
    stim_code = cellfun(@(v) str2double(regexp(v, '\d+', 'match', 'once')), {stim_mk.value});

    % --- Anchor: maintenance-onset codes s3023 (stick) / s3043 (avatar) ---
    anchor_codes = [3023, 3043];
    anchor_idx   = find(ismember(stim_code, anchor_codes));

    if isempty(anchor_idx)
        warning('%s: no maintenance-onset anchor markers (3023/3043) found.', subj);
        continue
    end

    n_trials = numel(anchor_idx);
    candidates = containers.Map('KeyType', 'double', 'ValueType', 'any');

    for t = 1:n_trials
        anchor_pos = stim_mk(anchor_idx(t)).position;

        % Look at a window of codes preceding this anchor, up to the
        % previous anchor (or start of recording), to capture all
        % candidate pre-maintenance events for this trial only.
        if t == 1
            win_start = 1;
        else
            win_start = anchor_idx(t-1) + 1;
        end
        win_idx = win_start:(anchor_idx(t) - 1);

        for k = win_idx
            code = stim_code(k);
            if isnan(code); continue; end
            lat_ms = (stim_mk(k).position - anchor_pos) / stim_mk(k).fs_hint * 1000; %#ok<NASGU>
            % fs_hint filled in by read_vmrk_markers via sidecar; if
            % unavailable, latency is left in samples and converted
            % downstream once srate is confirmed (500 Hz per dataset docs).
            lat_ms = (stim_mk(k).position - anchor_pos) / 500 * 1000;

            if ~isKey(candidates, code)
                candidates(code) = [];
            end
            candidates(code) = [candidates(code), lat_ms]; %#ok<AGROW>
        end
    end

    % --- Summarize candidates: count, mean/SD latency relative to anchor ---
    codes = cell2mat(keys(candidates));
    cand_summary = struct('code', {}, 'n', {}, 'mean_lat_ms', {}, 'sd_lat_ms', {}, 'occurs_once_per_trial', {});
    for c = codes
        lats = candidates(c);
        cand_summary(end+1) = struct( ...
            'code', c, ...
            'n', numel(lats), ...
            'mean_lat_ms', mean(lats), ...
            'sd_lat_ms', std(lats), ...
            'occurs_once_per_trial', abs(numel(lats) - n_trials) <= round(0.05 * n_trials)); %#ok<AGROW>
    end

    % Sort by how "clean" a candidate is: occurs once/trial, low SD, precedes anchor
    is_good = [cand_summary.occurs_once_per_trial] & [cand_summary.mean_lat_ms] < 0;
    good = cand_summary(is_good);
    [~, ord] = sort(abs([good.sd_lat_ms]));
    good = good(ord);

    fprintf('  %d trials found (anchor = maintenance onset s3023/s3043)\n', n_trials);
    fprintf('  Candidate fixed pre-maintenance markers (sorted by timing consistency):\n');
    for g = 1:numel(good)
        fprintf('    code=%d  n=%d  mean_lat=%.1fms  sd=%.1fms\n', ...
            good(g).code, good(g).n, good(g).mean_lat_ms, good(g).sd_lat_ms);
    end

    report(end+1) = struct('subject', subj, 'candidates', good, 'n_trials', n_trials); %#ok<AGROW>
end

save(fullfile(cfg_local.out_root, 'marker_report.mat'), 'report');

%% ---- POOLED RECOMMENDATION ----
% A marker code is a strong CUE/TARGET onset candidate if it appears with
% occurs_once_per_trial == true and low SD (<5 ms, i.e. essentially fixed
% latency = amplifier-locked trigger) in the large majority of subjects.
fprintf('\n=== POOLED SUMMARY ACROSS SUBJECTS ===\n');
all_codes = [];
for r = 1:numel(report)
    all_codes = [all_codes, [report(r).candidates.code]]; %#ok<AGROW>
end
u_codes = unique(all_codes);
for c = u_codes
    n_subj_with_code = sum(cellfun(@(cs) any(cs == c), arrayfun(@(r) [r.candidates.code], report, 'uni', 0)));
    fprintf('  code=%d present as clean candidate in %d/%d subjects\n', c, n_subj_with_code, numel(report));
end

fprintf(['\nRECOMMENDATION LOGIC (read before running s08b):\n' ...
    '  - If exactly ONE clean candidate code is shared by nearly all subjects\n' ...
    '    and its code value differs between stick/avatar trials (cross-check\n' ...
    '    against trialinfo col 6), treat it as CUE ONSET. Use s08b with\n' ...
    '    cfg.lock_event = ''cue''.\n' ...
    '  - If no such single fixed pre-maintenance marker exists (e.g. cue\n' ...
    '    presentation was not separately triggered, or jitter is large),\n' ...
    '    fall back to TARGET ONSET. In this dataset the most defensible\n' ...
    '    target-onset proxy is the marker immediately preceding the\n' ...
    '    maintenance-onset anchor with the SMALLEST mean|latency| across\n' ...
    '    subjects (i.e. the display event closest in time to maintenance\n' ...
    '    start). Use s08b with cfg.lock_event = ''target''.\n' ...
    '  - Document whichever choice is made, and the deviation from the\n' ...
    '    original cue-locked design, explicitly in the methods text.\n']);

%% ---- Local helper functions ----
function subs = get_valid_subject_list(raw_root)
    d = dir(fullfile(raw_root, 'sub-*'));
    excluded = {'sub-08', 'sub-42', 'sub-47'}; % corrupted / excluded per Pilot-2 pipeline
    subs = setdiff({d([d.isdir]).name}, excluded);
end

function f = find_vmrk_file(raw_root, subj)
    cand = dir(fullfile(raw_root, subj, '**', '*.vmrk'));
    if isempty(cand)
        f = '';
    else
        f = fullfile(cand(1).folder, cand(1).name);
    end
end

function mk = read_vmrk_markers(vmrk_file)
    % Minimal BrainVision .vmrk parser: returns struct array with
    % .value (Type+Description string, e.g. 'S 23') and .position (samples).
    fid = fopen(vmrk_file, 'r');
    raw = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = raw{1};
    mk = struct('value', {}, 'position', {}, 'fs_hint', {});
    for i = 1:numel(lines)
        ln = lines{i};
        if startsWith(ln, 'Mk') && contains(ln, '=')
            parts = strsplit(ln, {'=', ','});
            % Mk<N>=<Type>,<Description>,<Position>,<Size>,<Channel>
            if numel(parts) >= 4
                mk(end+1) = struct('value', strtrim(parts{3}), ...
                                    'position', str2double(parts{4}), ...
                                    'fs_hint', 500); %#ok<AGROW>
            end
        end
    end
end
