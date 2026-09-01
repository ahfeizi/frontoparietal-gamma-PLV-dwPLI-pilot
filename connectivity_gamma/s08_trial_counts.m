% s08_trial_counts.m  --  Trial-count table per subject x cell, post
% accuracy/RT exclusion (same filter as s02_epoch.m / s03_connectivity.m /
% s07_behavioral.m). Reports per-subject counts, aggregate stats, and flags
% subjects below a minimum-trial threshold -- relevant for gamma-band SNR
% and for dwPLI's (mild) sensitivity to unequal trial counts across cells.

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');   % has ntrials per subject/band/cue/valid

MIN_TRIALS = 15;   % reporting threshold; not a hard exclusion criterion

subs = unique(R.subject, 'stable');
cueLabs = {'avatar','stick'}; valLabs = {'valid','invalid'};

% ntrials is identical across bands (same trial-inclusion filter, band only
% affects the analytic signal) -- use gamma_low rows as the canonical count.
Rb = R(strcmp(R.band, 'gamma_low'), :);

rows = {};
for i = 1:numel(subs)
    subj = subs{i};
    row = {subj};
    for c = 1:numel(cueLabs)
        for v = 1:numel(valLabs)
            r = Rb(strcmp(Rb.subject,subj) & strcmp(Rb.cue,cueLabs{c}) & strcmp(Rb.valid,valLabs{v}), :);
            if isempty(r)
                row{end+1} = 0; %#ok<AGROW>
            else
                row{end+1} = r.ntrials(1); %#ok<AGROW>
            end
        end
    end
    rows(end+1,:) = row; %#ok<SAGROW>
end

T = cell2table(rows, 'VariableNames', {'subject','avatar_valid','avatar_invalid','stick_valid','stick_invalid'});
writetable(T, fullfile(cfg.deriv_dir, 'trial_counts_by_cell.csv'));

% ---- aggregate stats per cell ----
cellNames = {'avatar_valid','avatar_invalid','stick_valid','stick_invalid'};
fprintf('\n=== Trial counts per cell (N=%d subjects) ===\n', height(T));
fprintf('%-16s %8s %8s %8s %6s\n', 'cell', 'mean', 'min', 'max', sprintf('%%<%d', MIN_TRIALS));
summaryRows = {};
for i = 1:numel(cellNames)
    v = T.(cellNames{i});
    pctBelow = 100 * sum(v < MIN_TRIALS) / numel(v);
    fprintf('%-16s %8.1f %8d %8d %5.1f%%\n', cellNames{i}, mean(v), min(v), max(v), pctBelow);
    summaryRows(end+1,:) = {cellNames{i}, mean(v), std(v), min(v), max(v), pctBelow}; %#ok<SAGROW>
end
Summary = cell2table(summaryRows, 'VariableNames', {'cell','mean_n','sd_n','min_n','max_n','pct_below_thresh'});
writetable(Summary, fullfile(cfg.deriv_dir, 'trial_counts_summary.csv'));

% ---- overall (any-cell) flag ----
anyBelow = any(T{:,2:5} < MIN_TRIALS, 2);
fprintf('\nSubjects with at least one cell below %d trials: %d of %d (%.1f%%)\n', ...
    MIN_TRIALS, sum(anyBelow), height(T), 100*sum(anyBelow)/height(T));
if any(anyBelow)
    fprintf('Affected subjects:\n');
    disp(T(anyBelow,:));
end

% ---- valid vs invalid balance check (systematic imbalance flag) ----
diff_valid_invalid_avatar = T.avatar_valid - T.avatar_invalid;
diff_valid_invalid_stick  = T.stick_valid  - T.stick_invalid;
fprintf('\n=== Valid vs Invalid trial-count balance (paired sign-flip permutation) ===\n');
rng(1);
for lbl = {'avatar','stick'}
    if strcmp(lbl{1},'avatar'), d = diff_valid_invalid_avatar; else, d = diff_valid_invalid_stick; end
    d = d(~isnan(d)); n = numel(d);
    tobs = mean(d)/(std(d)/sqrt(n));
    np = cfg.stat.n_perm; tnull = zeros(np,1);
    for pp = 1:np
        dp = d .* sign(rand(n,1)-0.5);
        tnull(pp) = mean(dp)/(std(dp)/sqrt(n));
    end
    pval = mean(abs(tnull) >= abs(tobs));
    fprintf('%-8s: mean(valid-invalid) = %+.2f trials, dz=%.2f, t=%.2f, p=%.4f%s\n', ...
        lbl{1}, mean(d), mean(d)/std(d), tobs, pval, '');
    if pval < 0.05
        fprintf('  [FLAG] Significant trial-count imbalance between valid/invalid for %s -- \n', lbl{1});
        fprintf('         note as a potential dwPLI comparison caveat in Discussion.\n');
    end
end

fprintf('\nSaved trial_counts_by_cell.csv, trial_counts_summary.csv\n');
