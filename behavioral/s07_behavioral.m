% s07_behavioral.m  --  Behavioral check: accuracy and RT_loc as a function
% of cue (avatar/stick) x validity (valid/invalid), 2x2 within-subject,
% N=47. Uses the authors' trialinfo directly (cue=col6, validity=col4,
% accuracy=col16/17) with the SAME inclusion filter as s02_epoch.m
% (accuracy + RT_loc, with RT_stat disabled for the 43 subjects flagged in
% rt_qc_flags.csv).

clear; clc;
cfg = config();

sub_dirs = dir(fullfile(cfg.bids_dir, 'sub-*'));
sub_dirs = sub_dirs([sub_dirs.isdir]);

QC = readtable(fullfile(cfg.deriv_dir, 'rt_qc_flags.csv'));

rows = {};
for s = 1:numel(sub_dirs)
    subj = sub_dirs(s).name;
    ti_path = cfg.authors_deriv_file(subj);
    if ~exist(ti_path, 'file'), continue; end
    L = load(ti_path, 'trialinfo'); ti = L.trialinfo;

    cueV   = ti(:, cfg.ti.col_cue);
    validV = ti(:, cfg.ti.col_val);
    accV   = ti(:, cfg.ti.col_acc_loc);
    rt_loc = (ti(:, cfg.ti.col_locresp) - ti(:, cfg.ti.col_locprobe)) / cfg.resample_hz;

    qcRow = QC(strcmp(QC.subject, subj), :);
    statDegenerate = ~isempty(qcRow) && qcRow.rt_stat_degenerate(1);
    rt_stat = (ti(:, cfg.ti.col_statresp) - ti(:, cfg.ti.col_statq)) / cfg.resample_hz;

    rt_loc_ok = local_rt_ok(rt_loc, cfg);
    if statDegenerate
        rt_stat_ok = true(size(rt_loc));
    else
        rt_stat_ok = local_rt_ok(rt_stat, cfg);
    end

    good = (accV == cfg.accuracy_correct_value) & rt_loc_ok & rt_stat_ok;

    cueLabs = {'stick','avatar'}; cueValsArr = [cfg.ti.cue_stick, cfg.ti.cue_avatar];
    valLabs = {'valid','invalid'}; valValsArr = [cfg.ti.val_valid, cfg.ti.val_invalid];
    for c = 1:2
        cueLab = cueLabs{c}; cueVal = cueValsArr(c);
        for v = 1:2
            valLab = valLabs{v}; valVal = valValsArr(v);
            sel = good & cueV == cueVal & validV == valVal;
            acc_all = ti(cueV==cueVal & validV==valVal, cfg.ti.col_acc_loc);
            acc_cell = mean(acc_all, 'omitnan');
            rt_cell  = mean(rt_loc(sel), 'omitnan');
            rows(end+1,:) = {subj, cueLab, valLab, acc_cell, rt_cell, sum(sel)}; %#ok<SAGROW>
        end
    end
end

B = cell2table(rows, 'VariableNames', {'subject','cue','valid','accuracy','rt_loc_s','n_trials'});
writetable(B, fullfile(cfg.deriv_dir, 'behavioral_by_cell.csv'));

% ---- summary table: mean +- SD per cell ----
cueLabs = {'avatar','stick'}; valLabs = {'valid','invalid'};
fprintf('\n=== Summary: mean +/- SD per cell ===\n');
fprintf('%-16s %12s %14s\n', 'cell', 'accuracy', 'RT_loc (s)');
summaryRows = {};
for c = 1:2
    for v = 1:2
        sel = strcmp(B.cue,cueLabs{c}) & strcmp(B.valid,valLabs{v});
        acc = B.accuracy(sel); rt = B.rt_loc_s(sel);
        fprintf('%-16s %6.3f +/- %5.3f %6.3f +/- %5.3f\n', ...
            sprintf('%s-%s',cueLabs{c},valLabs{v}), mean(acc,'omitnan'), std(acc,'omitnan'), ...
            mean(rt,'omitnan'), std(rt,'omitnan'));
        summaryRows(end+1,:) = {sprintf('%s-%s',cueLabs{c},valLabs{v}), ...
            mean(acc,'omitnan'), std(acc,'omitnan'), mean(rt,'omitnan'), std(rt,'omitnan')}; %#ok<SAGROW>
    end
end
Summary = cell2table(summaryRows, 'VariableNames', {'cell','acc_mean','acc_sd','rt_mean','rt_sd'});
writetable(Summary, fullfile(cfg.deriv_dir, 'behavioral_summary.csv'));

% ---- 2x2 within-subject stats via sign-flip permutation (main effects + interaction) ----
rng(1);
subs = unique(B.subject, 'stable');
fprintf('\n=== 2x2 permutation tests ===\n');
fprintf('%-8s %-12s %3s %8s %8s %6s %6s %7s\n', 'DV','contrast','n','A','B','dz','t','p');
results = struct();
for dv = {'accuracy','rt_loc_s'}
    dvName = dv{1};
    av = local_cellval(B,subs,'avatar','valid',dvName);
    ai = local_cellval(B,subs,'avatar','invalid',dvName);
    sv = local_cellval(B,subs,'stick','valid',dvName);
    si = local_cellval(B,subs,'stick','invalid',dvName);

    contrasts = struct();
    contrasts.cue_main    = struct('A',(av+ai)/2,'B',(sv+si)/2);
    contrasts.valid_main  = struct('A',(av+sv)/2,'B',(ai+si)/2);
    contrasts.interaction = struct('A',av-ai,'B',sv-si);

    cnames = fieldnames(contrasts);
    for ci = 1:numel(cnames)
        cn = cnames{ci};
        A = contrasts.(cn).A; Bv = contrasts.(cn).B;
        ok = ~isnan(A) & ~isnan(Bv); d = A(ok)-Bv(ok); n = sum(ok);
        tobs = mean(d)/(std(d)/sqrt(n));
        np = cfg.stat.n_perm; tnull = zeros(np,1);
        for pp = 1:np
            dp = d .* sign(rand(n,1)-0.5);
            tnull(pp) = mean(dp)/(std(dp)/sqrt(n));
        end
        pval = mean(abs(tnull) >= abs(tobs));
        dzv = mean(d)/std(d);
        results.(dvName).(cn) = struct('n',n,'meanA',mean(A(ok)),'meanB',mean(Bv(ok)), ...
            'dz',dzv,'t',tobs,'p',pval);
        fprintf('%-8s %-12s %3d %8.4f %8.4f %6.2f %6.2f %.4f\n', ...
            dvName, cn, n, mean(A(ok)), mean(Bv(ok)), dzv, tobs, pval);
    end
end
save(fullfile(cfg.deriv_dir,'behavioral_stats.mat'), 'results', 'B', 'Summary');
fprintf('\nSaved behavioral_by_cell.csv, behavioral_summary.csv, behavioral_stats.mat\n');

% ======================= local functions =======================
function ok = local_rt_ok(rt, cfg)
    ok = ~isnan(rt) & (rt >= cfg.rt.min_rt_s);
    valid_rt = rt(ok);
    if numel(valid_rt) >= 5
        med = median(valid_rt); madv = median(abs(valid_rt - med));
        if madv > 0
            lo = med - cfg.rt.outlier_n*madv*1.4826;
            hi = med + cfg.rt.outlier_n*madv*1.4826;
            ok = ok & (rt >= lo) & (rt <= hi);
        end
    end
end

function v = local_cellval(B, subs, cueLab, valLab, dvName)
    n = numel(subs); v = nan(n,1);
    for i = 1:n
        r = B(strcmp(B.subject,subs{i}) & strcmp(B.cue,cueLab) & strcmp(B.valid,valLab), :);
        if ~isempty(r), v(i) = r.(dvName)(1); end
    end
end
