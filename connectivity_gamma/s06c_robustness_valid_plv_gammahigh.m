% s06c_robustness_valid_plv_gammahigh.m -- Robustness follow-up to
% s06_cluster_conn.m for the single finding it flagged: valid_main, PLV,
% gamma_high, pair F2-P4 (t=-3.696, p_fwer=0.0498, borderline significant
% under max-stat sign-flip permutation), sitting inside an apparently
% coherent frontal(right/mid)-parietal(right/mid) cluster of same-signed
% t-values that did not individually survive max-stat correction.
%
% Two checks, both restricted to valid_main / PLV / gamma_high to keep
% runtime bounded (can be generalized to other contrast/metric/band combos
% once the approach is confirmed useful):
%
%   (A) SEED STABILITY: s06 used a single fixed seed (rng(1)) and reused
%       cfg.stat.n_perm permutations. Because max-stat is a single
%       borderline p-value (0.0498), it is sensitive to permutation noise.
%       Re-run the identical max-stat test across several seeds, and again
%       with a larger n_perm, and report whether F2-P4 stays significant.
%
%   (B) CLUSTER-MASS TEST: max-stat (per-pair, Bonferroni-like via the max
%       null) is conservative for a spatially coherent effect and can miss
%       real clusters. This implements a standard cluster-mass permutation
%       test (Maris & Oostenveld, 2007): threshold each permutation's
%       per-pair t-map at a cluster-forming threshold, sum |t| within each
%       connected cluster, and compare the observed max cluster mass to the
%       null distribution of max cluster mass across permutations.
%
% *** ASSUMPTIONS THAT NEED CONFIRMATION -- NOT SILENTLY GUESSED ***
%   1. Cluster-forming threshold: not defined anywhere in cfg (s06 has no
%      such parameter). Using a parametric two-tailed p<0.05 threshold on
%      the observed t (df = n-1) purely to FORM clusters, which is standard
%      practice in Maris-Oostenveld even inside an otherwise nonparametric
%      pipeline -- but this specific alpha (0.05) is a new parameter and
%      should be confirmed/changed if a different convention is wanted.
%   2. Adjacency structure -- RESOLVED. Earlier versions used, in order:
%      (i) grid-index adjacency (invalidated by row-concatenated ROI lists),
%      (ii) NBS shared-endpoint adjacency (too coarse, no distance limit),
%      (iii) a heuristic 10-10-nomenclature-derived 2D layout (reasonable
%      but not verified against the actual pipeline's montage).
%      This version loads the REAL FieldTrip neighbour structure used by
%      the rest of the pipeline (neighbours.mat, built from layout_64WG.mat
%      via ft_prepare_neighbours) and restricts it to the frontal/parietal
%      ROI channels. Two frontal-parietal PAIRS are neighbors iff they
%      share an endpoint AND the other endpoint is a genuine montage
%      neighbor per neighbours.mat -- this is now authoritative, not
%      inferred.
%
% Path to the FieldTrip neighbour file.
clear; clc;
cfg = config();
neighbours_file = cfg.neighbours_file;
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');

nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_parietal); npair = nF*nT;
subs = unique(R.subject, 'stable');
cname = 'valid_main';
mt = 'plv';
bn = 'gamma_high';
fld = 'plv_pairs';

fprintf('=== Robustness check: %s / %s / %s ===\n', cname, mt, bn);

% ----- build the D matrix (subjects x pairs) for this band/metric/contrast -----
labs = cell(1,npair);
for i = 1:nF, for j = 1:nT
    labs{(i-1)*nT+j} = sprintf('%s-%s', cfg.roi_frontal{i}, cfg.roi_parietal{j});
end, end

D = nan(numel(subs), npair);
for si = 1:numel(subs)
    subj = subs{si};
    av = local_getpairs(R, subj, bn, 'avatar', 'valid',   fld);
    ai = local_getpairs(R, subj, bn, 'avatar', 'invalid', fld);
    sv = local_getpairs(R, subj, bn, 'stick',  'valid',   fld);
    si_ = local_getpairs(R, subj, bn, 'stick',  'invalid', fld);
    if isempty(av) || isempty(ai) || isempty(sv) || isempty(si_), continue; end
    D(si,:) = (av+sv)/2 - (ai+si_)/2;   % valid_main contrast
end
bad = any(isnan(D),2); D = D(~bad,:); n = size(D,1);
fprintf('Subjects used: %d\n\n', n);

tobs = mean(D) ./ (std(D)./sqrt(n));
target_idx = find(strcmp(labs, 'F2-P4'));
fprintf('Observed t for F2-P4: %.4f\n\n', tobs(target_idx));

% ======================= (A) seed stability, max-stat =======================
fprintf('--- (A) Max-stat seed stability ---\n');
seeds_to_try = [1 2 3 7 42];
n_perm_variants = [cfg.stat.n_perm, 10000];

fprintf('%8s %8s %10s %8s %8s\n','seed','n_perm','p_fwer(F2-P4)','sig','max|t|_thr');
for np = n_perm_variants
    for sd = seeds_to_try
        rng(sd);
        maxnull = zeros(np,1);
        for p = 1:np
            sgn = sign(rand(n,1)-0.5);
            Dp = D .* sgn;
            t = mean(Dp) ./ (std(Dp)./sqrt(n));
            maxnull(p) = max(abs(t));
        end
        p_fwer_target = mean(maxnull >= abs(tobs(target_idx)));
        thr95 = prctile(maxnull,95);
        fprintf('%8d %8d %10.4f %8d %8.3f\n', sd, np, p_fwer_target, p_fwer_target<0.05, thr95);
    end
end
fprintf(['\nInterpretation: if p_fwer(F2-P4) stays < 0.05 across seeds and\n' ...
    'grows more stable (converges) at n_perm=10000, the max-stat result is\n' ...
    'not a fluke of rng(1). If it flips above 0.05 for other seeds, the\n' ...
    'original p=0.0498 should be reported as non-robust / borderline.\n\n']);

% ======================= (B) cluster-mass permutation test =======================
fprintf('--- (B) Cluster-mass permutation test ---\n');
alpha_form = 0.05;                  % ASSUMPTION 1 -- confirm/change
v = n - 1;
t_thr = tinv(1 - alpha_form/2, v);  % two-tailed parametric cluster-forming threshold
fprintf('Cluster-forming threshold: |t| > %.3f (two-tailed p<%.2f, df=%d)\n', t_thr, alpha_form, v);

% ASSUMPTION 2 -- RESOLVED: real FieldTrip montage neighbours (see header
% note), loaded from neighbours.mat. Replaces the grid-index, NBS
% shared-endpoint, and 10-10-nomenclature-heuristic versions tried earlier.
adj = local_build_montage_adjacency(cfg.roi_frontal, cfg.roi_parietal, neighbours_file);

[obs_clusters, obs_mass] = local_find_clusters(tobs, adj, t_thr);
fprintf('Observed clusters at threshold: %d\n', numel(obs_mass));
for k = 1:numel(obs_clusters)
    fprintf('  cluster %d: %d pairs, mass=%.2f, members: %s\n', ...
        k, numel(obs_clusters{k}), obs_mass(k), strjoin(labs(obs_clusters{k}), ', '));
end

rng(1);
nperm = cfg.stat.n_perm;
max_mass_null = zeros(nperm,1);
for p = 1:nperm
    sgn = sign(rand(n,1)-0.5);
    Dp = D .* sgn;
    tp = mean(Dp) ./ (std(Dp)./sqrt(n));
    [~, mass_p] = local_find_clusters(tp, adj, t_thr);
    if isempty(mass_p), max_mass_null(p) = 0; else, max_mass_null(p) = max(mass_p); end
end

fprintf('\n%-10s %8s %10s %8s\n','cluster','mass','p_cluster','sig');
cluster_p = nan(numel(obs_mass),1);
for k = 1:numel(obs_mass)
    cluster_p(k) = mean(max_mass_null >= obs_mass(k));
    fprintf('%-10d %8.2f %10.4f %8d\n', k, obs_mass(k), cluster_p(k), cluster_p(k)<0.05);
end

T = table((1:numel(obs_mass))', obs_mass(:), cluster_p(:), cluster_p(:)<0.05, ...
    'VariableNames', {'cluster_id','mass','p_cluster','sig_cluster'});
writetable(T, fullfile(cfg.deriv_dir, sprintf('cluster_mass_%s_%s_%s.csv', cname, mt, bn)));
fprintf('\nSaved cluster_mass_%s_%s_%s.csv\n', cname, mt, bn);

% also save which pairs belong to which cluster, for the writeup
rows = {};
for k = 1:numel(obs_clusters)
    for idx = obs_clusters{k}
        rows(end+1,:) = {k, labs{idx}, tobs(idx)}; %#ok<SAGROW>
    end
end
Tmembers = cell2table(rows, 'VariableNames', {'cluster_id','pair','t'});
writetable(Tmembers, fullfile(cfg.deriv_dir, sprintf('cluster_mass_%s_%s_%s_members.csv', cname, mt, bn)));
fprintf('Saved cluster_mass_%s_%s_%s_members.csv\n', cname, mt, bn);

% ======================= local functions =======================
function pairs = local_getpairs(R, subj, bn, cueLab, valLab, fieldname)
    r = R(strcmp(R.subject,subj) & strcmp(R.band,bn) & ...
          strcmp(R.cue,cueLab) & strcmp(R.valid,valLab), :);
    if isempty(r)
        pairs = [];
    else
        pairs = r.(fieldname){1};
        if iscell(pairs), pairs = pairs{1}; end
    end
end

function adj = local_build_montage_adjacency(roi_frontal, roi_parietal, neighbours_file)
    S = load(neighbours_file, 'neighbours');
    nb = S.neighbours;
    nblabels = {nb.label};

    getneighbors = @(ch) local_lookup_neighbors(nb, nblabels, ch);

    nF = numel(roi_frontal); nT = numel(roi_parietal);
    npair = nF*nT;

    % restrict each channel's neighbor list to the ROI it belongs to
    Fadj = false(nF);
    for i = 1:nF
        nbrs = getneighbors(roi_frontal{i});
        for ip = 1:nF
            if ip~=i && any(strcmp(nbrs, roi_frontal{ip})), Fadj(i,ip) = true; end
        end
    end
    Padj = false(nT);
    for j = 1:nT
        nbrs = getneighbors(roi_parietal{j});
        for jp = 1:nT
            if jp~=j && any(strcmp(nbrs, roi_parietal{jp})), Padj(j,jp) = true; end
        end
    end

    adj = cell(npair,1);
    for i = 1:nF
        for j = 1:nT
            k = (i-1)*nT + j;
            nbrs_k = [];
            for jp = 1:nT
                if Padj(j,jp), nbrs_k(end+1) = (i-1)*nT + jp; end %#ok<AGROW>
            end
            for ip = 1:nF
                if Fadj(i,ip), nbrs_k(end+1) = (ip-1)*nT + j; end %#ok<AGROW>
            end
            adj{k} = unique(nbrs_k);
        end
    end
end

function nbrs = local_lookup_neighbors(nb, nblabels, ch)
    idx = find(strcmp(nblabels, ch), 1);
    if isempty(idx)
        warning('Channel %s not found in neighbours.mat -- treating as isolated (no neighbors). Verify ROI channel names match the montage exactly.', ch);
        nbrs = {};
    else
        nbrs = cellstr(nb(idx).neighblabel);
    end
end

function [clusters, mass] = local_find_clusters(tvec, adj, t_thr)
    npair = numel(tvec);
    supra = find(abs(tvec) > t_thr);
    visited = false(npair,1);
    clusters = {}; mass = [];
    for ii = 1:numel(supra)
        seed = supra(ii);
        if visited(seed), continue; end
        % BFS within same-signed supra-threshold pairs
        sgn = sign(tvec(seed));
        stack = seed; visited(seed) = true; comp = seed;
        while ~isempty(stack)
            cur = stack(end); stack(end) = [];
            for nb = adj{cur}
                if ~visited(nb) && abs(tvec(nb))>t_thr && sign(tvec(nb))==sgn
                    visited(nb) = true; stack(end+1) = nb; comp(end+1) = nb; %#ok<AGROW>
                end
            end
        end
        clusters{end+1} = comp; %#ok<AGROW>
        mass(end+1) = sum(abs(tvec(comp))); %#ok<AGROW>
    end
end
