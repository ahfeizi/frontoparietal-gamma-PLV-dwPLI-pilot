% s02b_check_maint_timing.m -- URGENT DIAGNOSTIC
% Checks whether the degenerate-RT_stat pattern found in 43/47 subjects
% (rt_qc_flags.csv) extends to the Maint->LocProbe interval, which is the
% ACTUAL connectivity analysis window (cfg.maint_window is defined relative
% to trialinfo column 11, Maint). If this interval is ALSO degenerate for
% the same subjects, our epoch definition itself is at risk and this
% becomes a much more serious issue than an unusable RT filter.

clear; clc;
cfg = config();

sub_dirs = dir(fullfile(cfg.bids_dir, 'sub-*'));
sub_dirs = sub_dirs([sub_dirs.isdir]);

rows = {};
for s = 1:numel(sub_dirs)
    subj = sub_dirs(s).name;
    ti_path = cfg.authors_deriv_file(subj);
    if ~exist(ti_path, 'file'), continue; end
    L = load(ti_path, 'trialinfo');
    ti = L.trialinfo;

    maint_to_probe = (ti(:, cfg.ti.col_locprobe) - ti(:, cfg.ti.col_maint)) / cfg.resample_hz;
    med = median(maint_to_probe, 'omitnan');
    mad_val = median(abs(maint_to_probe - med), 'omitnan');

    rows(end+1,:) = {subj, med, mad_val, min(maint_to_probe), max(maint_to_probe)}; %#ok<SAGROW>
    fprintf('%-10s Maint->LocProbe: median=%.4f s, MAD=%.4f s, range=[%.4f %.4f]\n', ...
        subj, med, mad_val, min(maint_to_probe), max(maint_to_probe));
end

T = cell2table(rows, 'VariableNames', {'subject','median_s','mad_s','min_s','max_s'});
writetable(T, fullfile(cfg.deriv_dir, 'maint_timing_check.csv'));

fprintf('\n%d subjects with MAD < 0.01 s (degenerate, same threshold as RT check): %d\n', ...
    sum(T.mad_s < cfg.rt.degenerate_mad_thresh_s), sum(T.mad_s < cfg.rt.degenerate_mad_thresh_s));
fprintf('Saved: %s\n', fullfile(cfg.deriv_dir, 'maint_timing_check.csv'));
