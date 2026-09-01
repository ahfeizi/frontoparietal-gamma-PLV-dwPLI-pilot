% s08_patch_add_elec.m  --  ONE-TIME PATCH, run only if you already ran
% s08_alpha_epoch_tfr.m BEFORE the elec-field fix and don't want to redo
% the (expensive) TFR computation for all subjects.
%
% Loads each subject's cleaned EEG just far enough to get electrode
% positions (via eeglab2fieldtrip), and appends subj_result.elec to the
% already-saved <key>_alphaTFR.mat file in place. Does NOT touch the TFR
% data itself. Idempotent: skips files that already have .elec.

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;

M = load(fullfile(cfg.deriv_dir, 'manifest.mat')); manifest = M.manifest;

n_patched = 0; n_skipped = 0; n_missing = 0;
for r = 1:height(manifest)
    key  = char(manifest.key(r));
    subj = char(manifest.subject(r));
    f = fullfile(cfg.deriv_dir, [key '_alphaTFR.mat']);

    if ~exist(f, 'file')
        n_missing = n_missing + 1;
        continue;
    end

    L = load(f);
    if isfield(L.subj_result, 'elec') && ~isempty(L.subj_result.elec)
        fprintf('%s already has elec, skipping.\n', key);
        n_skipped = n_skipped + 1;
        continue;
    end

    fprintf('Patching %s...\n', key);
    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    if isKey(cfg.manual_bad_channels, subj)
        EEG = pop_select(EEG, 'nochannel', cfg.manual_bad_channels(subj));
    end
    ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');

    subj_result = L.subj_result; %#ok<NASGU>
    subj_result.elec = ft_raw.elec;
    save(f, 'subj_result');
    n_patched = n_patched + 1;
end

fprintf('\nDone. Patched %d, already-ok %d, missing-file %d.\n', n_patched, n_skipped, n_missing);
fprintf('You can now run s09_alpha_stats.m.\n');
