% s16_multiple_comparisons_audit.m -- Multiple comparisons audit for
% Pilot 2, applying Holm-Bonferroni within two pre-specified families:
%
%   FAMILY 1 (behavioral, Phase 0): 4 DVs (acc_loc, acc_stat, rt_loc,
%   rt_stat) x 3 within-subject effects (cue main, valid main, cue:valid
%   interaction) = 12 tests. Source: behavior_anova_*.csv (s00b).
%
%   FAMILY 2 (S06 cluster FWER): 3 contrasts (cue_main, valid_main,
%   interaction) x 2 metrics (dwpli, plv) = 6 tests, each already
%   internally FWER-corrected across its own 128 pairs via max-stat
%   permutation (s06). Holm applied ACROSS these 6 already-corrected
%   family-level p-values (each family's own min p_fwer, i.e. its
%   strongest pair), since all 6 are reported together in the same
%   manuscript. Source: conn_cluster_results_<contrast>_by_metric.csv.
%
% EXCLUDED FROM FORMAL CORRECTION (by design, not oversight):
%   - S04 (24-test exploratory scan): superseded by S06/S10, not reported
%     as an independent statistical claim -- see manuscript note below.
%   - S10 (12-row Bayes/TOST table): Bayesian framework, no p-value alpha
%     to inflate; instead requires EXPLICIT DISCLOSURE that d_val/PLV/
%     gamma_high was singled out from this table for deeper follow-up
%     (S06c cluster-mass, S13 imaginary coherence) -- selection-after-
%     looking, not a pre-registered target.
%   - S06c / S13 (F2-P4 follow-up): confirmatory tests of an
%     already-identified target, not a new scan -- labeled "post-hoc" in
%     the manuscript rather than corrected as a new family.
%
% Output: multcomp_audit_family1_behavior.csv, multcomp_audit_family2_cluster.csv,
% plus manuscript-ready disclosure text printed to console.

clear; clc;
cfg = config();

% ======================= FAMILY 1: behavioral (Phase 0) =======================
fprintf('=== FAMILY 1: Phase 0 behavioral ANOVA (12 tests) ===\n');
dvs = {'acc_loc','acc_stat','rt_loc','rt_stat'};
effects = {'(Intercept):cue','(Intercept):valid','(Intercept):cue:valid'};
effect_labels = {'cue_main','valid_main','interaction'};

fam1_rows = {};
for d = 1:numel(dvs)
    fpath = fullfile(cfg.deriv_dir, sprintf('behavior_anova_%s.csv', dvs{d}));
    if ~isfile(fpath)
        warning('Missing %s -- skipping this DV in Family 1.', fpath);
        continue;
    end
    T = readtable(fpath);
    for e = 1:numel(effects)
        row = T(strcmp(T.term, effects{e}), :);
        if isempty(row)
            warning('Effect %s not found in %s -- skipping.', effects{e}, fpath);
            continue;
        end
        fam1_rows(end+1,:) = {dvs{d}, effect_labels{e}, row.F(1), row.pValue(1)}; %#ok<SAGROW>
    end
end
Fam1 = cell2table(fam1_rows, 'VariableNames', {'dv','effect','F','p_raw'});
Fam1 = local_holm(Fam1, 'p_raw');
writetable(Fam1, fullfile(cfg.deriv_dir, 'multcomp_audit_family1_behavior.csv'));
disp(Fam1);
fprintf('Saved multcomp_audit_family1_behavior.csv\n\n');

% ======================= FAMILY 2: S06 cluster (6 tests) =======================
fprintf('=== FAMILY 2: S06 cluster FWER, across-family Holm (6 tests) ===\n');
contrast_names = {'cue_main','valid_main','interaction'};
metrics = {'dwpli','plv'};
fam2_rows = {};
for c = 1:numel(contrast_names)
    fpath = fullfile(cfg.deriv_dir, sprintf('conn_cluster_results_%s_by_metric.csv', contrast_names{c}));
    if ~isfile(fpath)
        warning('Missing %s -- skipping this contrast in Family 2.', fpath);
        continue;
    end
    T = readtable(fpath);
    for m = 1:numel(metrics)
        Tm = T(strcmp(T.metric, metrics{m}), :);
        if isempty(Tm), continue; end
        % strongest (min) already-FWER-corrected p_fwer for this
        % contrast x metric combination -- the family-level result.
        [minp, idx] = min(Tm.p_fwer);
        fam2_rows(end+1,:) = {contrast_names{c}, metrics{m}, Tm.pair_band{idx}, Tm.t(idx), minp}; %#ok<SAGROW>
    end
end
Fam2 = cell2table(fam2_rows, 'VariableNames', {'contrast','metric','strongest_pair','t','p_fwer_raw'});
Fam2 = local_holm(Fam2, 'p_fwer_raw');
writetable(Fam2, fullfile(cfg.deriv_dir, 'multcomp_audit_family2_cluster.csv'));
disp(Fam2);
fprintf('Saved multcomp_audit_family2_cluster.csv\n\n');

% ======================= manuscript-ready disclosure text =======================
fprintf('=== Manuscript disclosure text (paste/adapt into Methods/Results) ===\n\n');
fprintf(['"An initial unc' 'orrected scan across band x metric x contrast combinations (s04-\n' ...
  'style) was used only to guide which analyses to pursue in depth and is not reported\n' ...
  'as an independent statistical claim. Two confirmatory families were subject to formal\n' ...
  'multiple-comparisons correction (Holm-Bonferroni): (1) the 12 behavioral ANOVA effects\n' ...
  '(4 DVs x 3 within-subject effects, Phase 0), and (2) the 6 cluster-level connectivity\n' ...
  'tests (3 contrasts x 2 metrics), each already FWER-corrected internally across its own\n' ...
  '128 channel pairs via max-statistic permutation, with Holm applied across the 6\n' ...
  'family-level results. The Bayesian gamma-band table (12 rows: 3 contrasts x 2 metrics x\n' ...
  '2 subbands) does not require frequentist alpha correction, but we explicitly disclose\n' ...
  'that the valid_main/PLV/gamma_high row was identified as notable from this table and\n' ...
  'selected for deeper follow-up (cluster-mass permutation with a montage-derived channel\n' ...
  'adjacency structure, and an imaginary-coherency check for volume-conduction artifacts);\n' ...
  'these follow-up analyses are confirmatory tests of an already-identified target, not a\n' ...
  'new scan, and are labeled post-hoc rather than corrected as an independent family."\n\n']);

% ======================= local functions =======================
function T = local_holm(T, pcol)
    [ps, order] = sort(T.(pcol));
    n = height(T);
    p_holm = nan(n,1);
    running_max = 0;
    for i = 1:n
        adj = (n - i + 1) * ps(i);
        running_max = max(running_max, adj);
        p_holm(order(i)) = min(running_max, 1);
    end
    T.p_holm = p_holm;
    T.sig_holm_05 = p_holm < 0.05;
end
