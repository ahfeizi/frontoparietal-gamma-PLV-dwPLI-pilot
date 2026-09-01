% s00b_phase0_behavior.m -- Phase 0 for Pilot 2: behavioral ANOVA and
% trial-count/retention table, BEFORE any EEG connectivity analysis.
% Built directly from each subject's <key>_trials.mat (s02_epoch.m output),
% NOT from raw trialinfo, so it reflects exactly what the EEG analyses use.
%
% Verified before writing: acc_loc and acc_stat in trials.mat are confirmed
% binary {0,1} (0=incorrect, 1=correct) on real data (sub-01), matching
% config.m's documented convention -- not just trusted from the config
% comment, but re-checked on the actual file the rest of the pipeline reads.
%
% Outputs:
%   1. trial_count_table.csv   -- per subject x cue x validity: n_raw,
%      n_good (post-QC), pct_retained
%   2. behavior_anova_acc_loc.csv / behavior_anova_acc_stat.csv -- 2x2
%      within-subject ANOVA (cue x validity) on location/status accuracy
%   3. behavior_rt_summary.csv -- descriptive + ANOVA on rt_loc and rt_stat,
%      with an explicit caveat column flagging subjects whose RT_stat was
%      auto-disabled by the QC pipeline (rt_qc_flags.csv, if present)

clear; clc;
cfg = config();
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

cue_code   = struct('avatar', cfg.ti.cue_avatar, 'stick', cfg.ti.cue_stick);
valid_code = struct('valid',  cfg.ti.val_valid,  'invalid', cfg.ti.val_invalid);
cueLabs = {'avatar','stick'}; valLabs = {'valid','invalid'};

% Optional: RT_stat QC flags from the connectivity pipeline (per your
% earlier notes: 43/47 subjects had degenerate RT_stat, auto-disabled).
% Loaded if present so the RT summary can flag it explicitly rather than
% silently mixing flagged and unflagged subjects.
rt_flag_file = fullfile(cfg.deriv_dir, 'rt_qc_flags.csv');
have_rt_flags = isfile(rt_flag_file);
if have_rt_flags
    RTflags = readtable(rt_flag_file);
else
    warning('rt_qc_flags.csv not found at %s -- RT summary will not mark degenerate-RT subjects.', rt_flag_file);
end

count_rows = {};
acc_loc_rows = {}; acc_stat_rows = {};
rt_loc_rows = {}; rt_stat_rows = {};

for su = 1:height(manifest)
    subj = char(manifest.subject(su));
    key  = char(manifest.key(su));

    trials_file = fullfile(cfg.deriv_dir, [key '_trials.mat']);
    if ~isfile(trials_file)
        warning('%s: trials file not found (%s), skipping.', subj, trials_file);
        continue;
    end
    T = load(trials_file); trials = T.trials;

    cueArr = [trials.cue]; validArr = [trials.valid]; goodArr = logical([trials.good]);
    accLocArr = [trials.acc_loc]; accStatArr = [trials.acc_stat];
    rtLocArr = [trials.rt_loc]; rtStatArr = [trials.rt_stat];

    for c = 1:numel(cueLabs)
        for v = 1:numel(valLabs)
            cueLab = cueLabs{c}; valLab = valLabs{v};
            sel_raw = (cueArr==cue_code.(cueLab)) & (validArr==valid_code.(valLab));
            sel_good = sel_raw & goodArr;

            n_raw = sum(sel_raw); n_good = sum(sel_good);
            pct = 100 * n_good / max(n_raw,1);
            count_rows(end+1,:) = {subj, cueLab, valLab, n_raw, n_good, pct}; %#ok<SAGROW>

            % accuracy: computed on ALL RAW trials for this cell, NOT on
            % "good" trials. cfg.correct_only=true means the accuracy
            % criterion is already baked into the "good" flag (per
            % config.m), so restricting to good trials would make acc_loc
            % trivially 1.0 for every subject/condition (confirmed: this
            % is exactly what happened on the first run). Using sel_raw
            % matches the original paper's approach (accuracy computed
            % over all attempted trials, not a QC-filtered subset).
            acc_loc_rows(end+1,:)  = {subj, cueLab, valLab, mean(accLocArr(sel_raw), 'omitnan')}; %#ok<SAGROW>
            acc_stat_rows(end+1,:) = {subj, cueLab, valLab, mean(accStatArr(sel_raw), 'omitnan')}; %#ok<SAGROW>

            % RT: computed on GOOD trials (i.e. after RT-outlier/anticipatory
            % exclusion), which is standard practice (RT on correct,
            % non-outlier trials only) -- NOTE this means the RT distribution
            % here is truncated relative to a fully unfiltered RT, which is
            % intentional but should be stated in the manuscript.
            rt_loc_rows(end+1,:)   = {subj, cueLab, valLab, mean(rtLocArr(sel_good), 'omitnan'), sum(~isnan(rtLocArr(sel_good)))}; %#ok<SAGROW>

            rt_stat_flagged = false;
            if have_rt_flags && any(strcmp(RTflags.subject, subj))
                rowf = RTflags(strcmp(RTflags.subject,subj),:);
                if any(rowf.rt_stat_degenerate == 1)   % CONFIRMED column name/values
                    rt_stat_flagged = true;             % from rt_qc_flags.csv schema
                end
            end
            rt_stat_rows(end+1,:) = {subj, cueLab, valLab, mean(rtStatArr(sel_good), 'omitnan'), ...
                sum(~isnan(rtStatArr(sel_good))), rt_stat_flagged}; %#ok<SAGROW>
        end
    end
end

Tcount = cell2table(count_rows, 'VariableNames', {'subject','cue','valid','n_raw','n_good','pct_retained'});
writetable(Tcount, fullfile(cfg.deriv_dir,'trial_count_table.csv'));
fprintf('Saved trial_count_table.csv (%d subject x cell rows)\n', height(Tcount));


% ======================= trial-count summary printed to console =======================
fprintf('\n--- Trial retention summary (across subjects x cells) ---\n');
fprintf('n_raw:  median=%.1f, range=[%d %d]\n', median(Tcount.n_raw), min(Tcount.n_raw), max(Tcount.n_raw));
fprintf('n_good: median=%.1f, range=[%d %d]\n', median(Tcount.n_good), min(Tcount.n_good), max(Tcount.n_good));
fprintf('pct_retained: median=%.1f%%, range=[%.1f%% %.1f%%]\n', ...
    median(Tcount.pct_retained), min(Tcount.pct_retained), max(Tcount.pct_retained));
low_cells = Tcount(Tcount.n_good < 5, :);
if ~isempty(low_cells)
    fprintf('\nWARNING: %d subject x cell combinations have < 5 good trials:\n', height(low_cells));
    disp(low_cells);
else
    fprintf('No subject x cell combination has fewer than 5 good trials.\n');
end

fprintf(['\nNOTE: acc_loc/acc_stat below are computed on ALL raw trials per cell ' ...
    '(not the QC-filtered "good" set), because cfg.correct_only=true bakes the ' ...
    'accuracy criterion into "good" -- computing accuracy on "good" trials only ' ...
    'would be circular (confirmed empirically: acc_loc=1.0 for every subject/\n' ...
    'condition on the first run). rt_loc/rt_stat remain computed on "good" trials\n' ...
    '(RT-outlier-filtered), matching standard practice of analyzing RT on valid\n' ...
    'responses only -- this truncation should be noted in the manuscript.\n\n']);

% ======================= behavioral ANOVA: acc_loc =======================
Tal = cell2table(acc_loc_rows, 'VariableNames', {'subject','cue','valid','acc'});
run_2x2_anova(Tal, 'acc_loc', cfg.deriv_dir);

% ======================= behavioral ANOVA: acc_stat =======================
Tas = cell2table(acc_stat_rows, 'VariableNames', {'subject','cue','valid','acc'});
run_2x2_anova(Tas, 'acc_stat', cfg.deriv_dir);

% ======================= RT summary + ANOVA (both response types) =======================
Trl = cell2table(rt_loc_rows, 'VariableNames', {'subject','cue','valid','rt','n_valid_rt'});
run_2x2_anova(Trl, 'rt_loc', cfg.deriv_dir);

Trs = cell2table(rt_stat_rows, 'VariableNames', {'subject','cue','valid','rt','n_valid_rt','rt_stat_flagged'});
writetable(Trs, fullfile(cfg.deriv_dir,'behavior_rt_stat_raw.csv'));
n_flagged_subj = numel(unique(Trs.subject(Trs.rt_stat_flagged)));
fprintf(['\nCAVEAT: rt_stat ANOVA/summary includes %d subjects flagged elsewhere as having ' ...
    'degenerate RT_stat (near-constant values, logged in rt_qc_flags.csv). Their rt_stat ' ...
    'means are included below for completeness but should be reported/interpreted with this ' ...
    'caveat explicit in the manuscript -- do not present rt_stat effects without noting it.\n'], n_flagged_subj);
run_2x2_anova(Trs(:,{'subject','cue','valid','rt'}), 'rt_stat', cfg.deriv_dir);

% ======================= local functions =======================
function run_2x2_anova(Tlong, label, deriv_dir)
    % Reshapes a long (subject,cue,valid,DV) table into a 2x2 within-subject
    % design and runs a repeated-measures ANOVA via MATLAB's fitrm/ranova.
    subs = unique(Tlong.subject, 'stable');
    wide = table(subs, 'VariableNames', {'subject'});
    condnames = {'avatar_valid','avatar_invalid','stick_valid','stick_invalid'};
    combos = {'avatar','valid'; 'avatar','invalid'; 'stick','valid'; 'stick','invalid'};
    for k = 1:4
        v = nan(numel(subs),1);
        for i = 1:numel(subs)
            r = Tlong(strcmp(Tlong.subject,subs{i}) & strcmp(Tlong.cue,combos{k,1}) & strcmp(Tlong.valid,combos{k,2}), :);
            if ~isempty(r), v(i) = r.(Tlong.Properties.VariableNames{4})(1); end
        end
        wide.(condnames{k}) = v;
    end
    bad = any(isnan(wide{:,2:5}),2);
    nDrop = sum(bad);
    wide = wide(~bad,:);
    if nDrop > 0
        fprintf('%s ANOVA: dropped %d/%d subjects with missing cells.\n', label, nDrop, numel(subs));
    end

    % DIAGNOSTIC: print raw cell values for the first few subjects before
    % running the ANOVA, so a degenerate/near-constant DV (e.g. ceiling
    % accuracy, or a data-extraction bug) is visible rather than only
    % showing up as a suspicious all-zero-variance ANOVA table.
    fprintf('\n--- %s: raw per-subject condition means (first 8 subjects) ---\n', label);
    disp(wide(1:min(8,height(wide)),:));
    fprintf('%s overall: mean=%.4f, sd across subjects (grand-mean per subject)=%.6f\n', ...
        label, mean(wide{:,2:5}(:)), std(mean(wide{:,2:5},2)));

    within = table(categorical({'avatar';'avatar';'stick';'stick'}), ...
                   categorical({'valid';'invalid';'valid';'invalid'}), ...
                   'VariableNames', {'cue','valid'});
    rm = fitrm(wide, sprintf('%s-%s ~ 1', condnames{1}, condnames{4}), 'WithinDesign', within);
    ranovatbl = ranova(rm, 'WithinModel', 'cue*valid');
    disp(label); disp(ranovatbl);

    % ranova() returns a special anova-table class whose row names are
    % compound labels (e.g. "(Intercept):cue"); writetable(...,'WriteRowNames',true)
    % on this class throws an internal error ("DoubleTableColumnReal") in
    % some MATLAB versions. Convert to a plain table with an explicit
    % 'term' column instead, which is both safer and clearer to read.
    plainT = table(string(ranovatbl.Properties.RowNames), 'VariableNames', {'term'});
    varnames = ranovatbl.Properties.VariableNames;
    for vn = 1:numel(varnames)
        col = ranovatbl.(varnames{vn});
        try
            col = double(col);   % strip any internal display-wrapper class (e.g.
        catch                    % internal.stats.DoubleTableColumnReal) that can't
        end                      % be assigned directly into a plain table column
        plainT.(varnames{vn}) = col;
    end
    writetable(plainT, fullfile(deriv_dir, sprintf('behavior_anova_%s.csv', label)));

    % descriptive means per cell, for the manuscript table
    means = mean(wide{:,2:5}, 1); sds = std(wide{:,2:5}, 0, 1);
    fprintf('  cell means (%s): avatar_valid=%.3f(%.3f) avatar_invalid=%.3f(%.3f) stick_valid=%.3f(%.3f) stick_invalid=%.3f(%.3f)\n', ...
        label, means(1),sds(1), means(2),sds(2), means(3),sds(3), means(4),sds(4));
end
