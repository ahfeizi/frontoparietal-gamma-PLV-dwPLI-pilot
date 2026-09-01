% s20_csd_imc_cluster_followup.m -- Two follow-up checks on the new
% 9-pair ImC/CSD cluster found by s19_csd_perpair_imc.m for valid_main/
% gamma_high (mass=24.41, p_cluster=0.0036, members: F1-P3, Fz-P3, Fz-P1,
% Fz-Pz, Fz-POz, F1-POz, Fz-P2, FCz-POz, FC1-P2):
%
%   (A) SEED STABILITY of the cluster-mass permutation test. CSD +
%       connectivity is recomputed ONCE (the expensive part, ~1-2 min for
%       47 subjects, identical recipe to s19) and the resulting per-
%       subject/per-pair D matrix (d_val contrast, ImC) is CACHED to a
%       .mat file. The permutation test itself (cheap) is then re-run
%       across 5 seeds x 2 n_perm values, tracking p_cluster for THIS
%       specific cluster's observed mass (24.41) each time -- same
%       stability-check logic used earlier for the F2-P4 max-stat result
%       (s06c), adapted to cluster-mass.
%
%   (B) CROSS-METRIC CHECK: are dwPLI/PLV (also computed and saved by
%       s19, in csd_imc_perpair_valid_main_gamma_high.csv) elevated -- even
%       if not significant -- at the same 9 pairs? This does NOT require
%       new computation; it just filters the already-saved per-pair table.
%
% Output:
%   csd_imc_cluster_seedstability.csv
%   csd_imc_cluster_crossmetric.csv

clear; clc;
cfg = config();

target_pairs = {'F1-P3','Fz-P3','Fz-P1','Fz-POz','F1-POz','Fz-P2','FCz-POz'};
observed_mass = 17.12;   % from s19's console output (post s01 muscle-IC fix), for
                          % reference/sanity-check only; the script recomputes its own
                          % observed mass from the cached D matrix rather than trusting
                          % this hardcoded value for the actual test.

% ======================= (B) cross-metric check -- cheap, do first ======================
fprintf('=== (B) Cross-metric check: dwPLI/PLV at the 9 ImC-cluster pairs (CSD-referenced) ===\n');
perpair_file = fullfile(cfg.deriv_dir, 'csd_imc_perpair_valid_main_gamma_high.csv');
if ~isfile(perpair_file)
    warning('%s not found -- run s19_csd_perpair_imc.m first. Skipping (B).', perpair_file);
else
    Tpp = readtable(perpair_file);
    Tsub = Tpp(ismember(Tpp.pair, target_pairs) & ismember(Tpp.metric, {'dwPLI','PLV'}), :);
    Tsub = sortrows(Tsub, {'metric','pair'});
    disp(Tsub(:, {'pair','metric','n','dz','t_freq','p_freq','BF10','BF01'}));
    writetable(Tsub, fullfile(cfg.deriv_dir, 'csd_imc_cluster_crossmetric.csv'));
    fprintf('Saved csd_imc_cluster_crossmetric.csv\n');
    n_dwpli_nominal = sum(strcmp(Tsub.metric,'dwPLI') & Tsub.p_freq<0.05);
    n_plv_nominal = sum(strcmp(Tsub.metric,'PLV') & Tsub.p_freq<0.05);
    fprintf('Of 9 pairs: %d/9 nominally significant (p<.05) in dwPLI, %d/9 in PLV.\n', ...
        n_dwpli_nominal, n_plv_nominal);
end

% ======================= (A) seed stability of the cluster-mass test =======================
fprintf('\n=== (A) Seed stability: recomputing CSD + ImC once, caching D matrix ===\n');
cache_file = fullfile(cfg.deriv_dir, 'csd_imc_dmatrix_cache.mat');

% Invalidate the cache if it's older than the most recently modified
% _clean.set file, so a preprocessing rerun can never silently produce
% stale downstream numbers again (this is exactly the bug that made this
% script report mass=24.41/9 pairs from the pre-fix data even after s19
% had already been correctly rerun on post-fix data).
if isfile(cache_file)
    cache_info = dir(cache_file);
    set_files = dir(fullfile(cfg.deriv_dir, '*_clean.set'));
    newest_set = max([set_files.datenum]);
    if cache_info.datenum < newest_set
        warning('csd_imc_dmatrix_cache.mat is older than the newest _clean.set file -- deleting stale cache.');
        delete(cache_file);
    end
end

if isfile(cache_file)
    fprintf('Found cached D matrix (%s) -- reusing instead of recomputing CSD.\n', cache_file);
    L = load(cache_file); D = L.D; labs = L.labs; subj_list = L.subj_list;
else
    addpath(cfg.eeglab_dir); eeglab nogui;
    addpath(cfg.fieldtrip_dir); ft_defaults;
    M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

    bn = 'gamma_high'; band_range = cfg.bands.(bn);
    srate = cfg.resample_hz; win = round(cfg.maint_window * srate); wlen = win(2)-win(1)+1;
    nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_parietal);
    labs = cell(1,nF*nT);
    for i=1:nF, for j=1:nT, labs{(i-1)*nT+j} = sprintf('%s-%s', cfg.roi_frontal{i}, cfg.roi_parietal{j}); end, end
    cueVals = [cfg.ti.cue_avatar, cfg.ti.cue_stick]; valVals = [cfg.ti.val_valid, cfg.ti.val_invalid];

    nSub = height(manifest);
    IMC = nan(nSub, 4, nF*nT); subj_list = cell(nSub,1);

    for su = 1:nSub
        subj = char(manifest.subject(su)); key = char(manifest.key(su));
        subj_list{su} = subj;
        trials_file = fullfile(cfg.deriv_dir, [key '_trials.mat']);
        set_file = fullfile(cfg.deriv_dir, [key '_clean.set']);
        if ~isfile(trials_file) || ~isfile(set_file), continue; end
        T = load(trials_file); trials = T.trials;
        good_trials = trials([trials.good]);
        if isempty(good_trials), continue; end

        EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
        ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
        cfg_csd = []; cfg_csd.method = 'spline';
        ft_csd = ft_scalpcurrentdensity(cfg_csd, ft_raw);
        csd_data = ft_csd.trial{1}; csd_labels = ft_csd.label;

        fr = local_chanidx(csd_labels, cfg.roi_frontal);
        pa = local_chanidx(csd_labels, cfg.roi_parietal);
        if numel(fr) < numel(cfg.roi_frontal) || numel(pa) < numel(cfg.roi_parietal), continue; end

        Xf = util_bandpass(double(csd_data), srate, band_range, cfg.filt_order_factor);
        Z  = hilbert(Xf.').';

        for c = 1:2
            for v = 1:2
                cidx = (c-1)*2 + v;
                sel = good_trials([good_trials.cue]==cueVals(c) & [good_trials.valid]==valVals(v));
                ntr = numel(sel);
                if ntr < 5, continue; end
                seg = nan(numel(csd_labels), wlen, ntr); kept = 0;
                for t = 1:ntr
                    a = sel(t).lock_sample + win(1); z = a + wlen - 1;
                    if a < 1 || z > size(Z,2), continue; end
                    kept = kept+1; seg(:,:,kept) = Z(:, a:z);
                end
                seg = seg(:,:,1:kept);
                if kept < 5, continue; end
                p = 0;
                for i = 1:numel(fr)
                    zi = reshape(seg(fr(i),:,:), wlen, kept);
                    for j = 1:numel(pa)
                        p = p+1;
                        zj = reshape(seg(pa(j),:,:), wlen, kept);
                        X = zi.*conj(zj);
                        Sxy = mean(mean(X)); Sxx = mean(mean(abs(zi).^2)); Syy = mean(mean(abs(zj).^2));
                        IMC(su,cidx,p) = imag(Sxy/sqrt(Sxx*Syy));
                    end
                end
            end
        end
        fprintf('  %s done\n', subj);
    end

    av = IMC(:,1,:); ai = IMC(:,2,:); sv = IMC(:,3,:); iv = IMC(:,4,:);
    D = squeeze((av+sv)/2 - (ai+iv)/2);
    save(cache_file, 'D', 'labs', 'subj_list');
    fprintf('Saved %s for reuse in future seed/threshold checks.\n', cache_file);
end

bad = any(isnan(D),2); D = D(~bad,:); n = size(D,1);
tobs = mean(D) ./ (std(D)./sqrt(n));
target_idx = find(ismember(labs, target_pairs));
fprintf('n=%d subjects.\n', n);
neighbours_file = cfg.neighbours_file;
adj = local_build_montage_adjacency(cfg.roi_frontal, cfg.roi_parietal, neighbours_file);
alpha_form = 0.05; v = n-1; t_thr = tinv(1-alpha_form/2, v);

[obs_clusters, obs_mass] = local_find_clusters(tobs, adj, t_thr);
k_target = find(cellfun(@(c) numel(intersect(c,target_idx))>=numel(target_idx)/2, obs_clusters), 1);
if isempty(k_target)
    error('Could not re-identify the target 9-pair cluster in the recomputed D matrix -- membership may have shifted. Inspect obs_clusters manually before trusting seed-stability numbers.');
end
fprintf('Re-identified target cluster: %d pairs, observed mass=%.2f (s19 reported %.2f)\n', ...
    numel(obs_clusters{k_target}), obs_mass(k_target), observed_mass);

seeds_to_try = [1 2 3 7 42];
n_perm_variants = [cfg.stat.n_perm, 10000];
results_rows = {};
fprintf('\n%8s %8s %10s %8s\n','seed','n_perm','mass','p_cluster');
for np = n_perm_variants
    for sd = seeds_to_try
        rng(sd);
        max_mass_null = zeros(np,1);
        for pp = 1:np
            sgn = sign(rand(n,1)-0.5); Dp = D.*sgn; tp = mean(Dp)./(std(Dp)./sqrt(n));
            [~, mass_p] = local_find_clusters(tp, adj, t_thr);
            if isempty(mass_p), max_mass_null(pp)=0; else, max_mass_null(pp)=max(mass_p); end
        end
        p_clust = mean(max_mass_null >= obs_mass(k_target));
        fprintf('%8d %8d %10.4f %8.4f\n', sd, np, obs_mass(k_target), p_clust);
        results_rows(end+1,:) = {sd, np, obs_mass(k_target), p_clust}; %#ok<SAGROW>
    end
end
Tseed = cell2table(results_rows, 'VariableNames', {'seed','n_perm','mass','p_cluster'});
writetable(Tseed, fullfile(cfg.deriv_dir, 'csd_imc_cluster_seedstability.csv'));
fprintf('\nSaved csd_imc_cluster_seedstability.csv\n');
fprintf(['\nInterpretation: if p_cluster stays comfortably below 0.05 across all seeds and\n' ...
    'both n_perm values, the cluster is robust and not an artifact of rng(1). If it is\n' ...
    'unstable (crossing 0.05), report the original 0.0036 as non-robust / requiring a\n' ...
    'larger n_perm as the default going forward.\n']);

% ======================= local functions =======================
function idx = local_chanidx(labels, want)
    idx = [];
    for i = 1:numel(want)
        j = find(strcmpi(labels, want{i}), 1);
        if ~isempty(j), idx(end+1) = j; end %#ok<AGROW>
    end
end

function adj = local_build_montage_adjacency(roi_frontal, roi_parietal, neighbours_file)
    S = load(neighbours_file, 'neighbours'); nb = S.neighbours; nblabels = {nb.label};
    getneighbors = @(ch) local_lookup(nb, nblabels, ch);
    nF=numel(roi_frontal); nT=numel(roi_parietal); npair=nF*nT;
    Fadj=false(nF); for i=1:nF, nbrs=getneighbors(roi_frontal{i}); for ip=1:nF, if ip~=i && any(strcmp(nbrs,roi_frontal{ip})), Fadj(i,ip)=true; end, end, end
    Padj=false(nT); for j=1:nT, nbrs=getneighbors(roi_parietal{j}); for jp=1:nT, if jp~=j && any(strcmp(nbrs,roi_parietal{jp})), Padj(j,jp)=true; end, end, end
    adj=cell(npair,1);
    for i=1:nF, for j=1:nT
        k=(i-1)*nT+j; nbrs_k=[];
        for jp=1:nT, if Padj(j,jp), nbrs_k(end+1)=(i-1)*nT+jp; end, end %#ok<AGROW>
        for ip=1:nF, if Fadj(i,ip), nbrs_k(end+1)=(ip-1)*nT+j; end, end %#ok<AGROW>
        adj{k}=unique(nbrs_k);
    end, end
end

function nbrs = local_lookup(nb, nblabels, ch)
    idx = find(strcmp(nblabels,ch),1);
    if isempty(idx), nbrs={}; else, nbrs=cellstr(nb(idx).neighblabel); end
end

function [clusters, mass] = local_find_clusters(tvec, adj, t_thr)
    npair=numel(tvec); supra=find(abs(tvec)>t_thr); visited=false(npair,1);
    clusters={}; mass=[];
    for ii=1:numel(supra)
        seed=supra(ii); if visited(seed), continue; end
        sgn=sign(tvec(seed)); stack=seed; visited(seed)=true; comp=seed;
        while ~isempty(stack)
            cur=stack(end); stack(end)=[];
            for nb=adj{cur}
                if ~visited(nb) && abs(tvec(nb))>t_thr && sign(tvec(nb))==sgn
                    visited(nb)=true; stack(end+1)=nb; comp(end+1)=nb; %#ok<AGROW>
                end
            end
        end
        clusters{end+1}=comp; mass(end+1)=sum(abs(tvec(comp))); %#ok<AGROW>
    end
end
