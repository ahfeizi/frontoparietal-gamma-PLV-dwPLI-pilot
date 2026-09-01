% s19_csd_perpair_imc.m -- Follow-up to s18_csd_sensitivity.m: the
% CSD-referenced ImC ROI-average result for valid_main/gamma_high showed
% a moderate, unexpected effect (BF10=3.74, p_freq=0.010) not present
% under average reference (where ImC was null: BF01=6.18) and not present
% at the single F2-P4 pair under CSD either (BF01=6.08). This script
% determines whether that ROI-average anomaly is localized to a specific
% subset of frontal-parietal pairs (consistent with a real, CSD-revealed
% weak signal) or diffuse/uniform across all 64 pairs (more consistent
% with CSD amplifying sensor noise that happens to average out non-zero
% by chance, or a multiple-comparisons artifact).
%
% Recomputes CSD + Imaginary Coherency for valid_main/gamma_high, per
% pair, for all 47 subjects (same recipe as s18: branch from the already
% saved <key>_clean.set via ft_scalpcurrentdensity('spline'), justified by
% the reference-invariance of the surface Laplacian -- see s18 header).
% dwPLI/PLV are also recomputed in the same loop (negligible extra cost)
% and included per-pair for completeness/consistency checking, though the
% anomaly under investigation is specifically in ImC.
%
% Output:
%   csd_imc_perpair_valid_main_gamma_high.csv -- 64 rows, BF/TOST/freq
%     stats per pair, for dwPLI, PLV, and ImC (192 rows total: 64 x 3)
%   csd_imc_clustermass_valid_main_gamma_high.csv -- cluster-mass results
%     for ImC specifically, using the real montage adjacency (neighbours.mat)

clear; clc;
cfg = config();
addpath(cfg.eeglab_dir); eeglab nogui;
addpath(cfg.fieldtrip_dir); ft_defaults;

M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

bn = 'gamma_high';
band_range = cfg.bands.(bn);
srate = cfg.resample_hz;
win = round(cfg.maint_window * srate);
wlen = win(2) - win(1) + 1;

nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_parietal);
labs = cell(1,nF*nT);
for i=1:nF, for j=1:nT, labs{(i-1)*nT+j} = sprintf('%s-%s', cfg.roi_frontal{i}, cfg.roi_parietal{j}); end, end

cueVals = [cfg.ti.cue_avatar, cfg.ti.cue_stick];
valVals = [cfg.ti.val_valid, cfg.ti.val_invalid];

nSub = height(manifest);
DWPLI = nan(nSub, 4, nF*nT); PLV = nan(nSub, 4, nF*nT); IMC = nan(nSub, 4, nF*nT);
subj_list = cell(nSub,1);

for su = 1:nSub
    subj = char(manifest.subject(su)); key = char(manifest.key(su));
    subj_list{su} = subj;

    trials_file = fullfile(cfg.deriv_dir, [key '_trials.mat']);
    set_file = fullfile(cfg.deriv_dir, [key '_clean.set']);
    if ~isfile(trials_file) || ~isfile(set_file)
        warning('%s: missing trials/clean.set, skipping.', subj); continue;
    end
    T = load(trials_file); trials = T.trials;
    good_trials = trials([trials.good]);
    if isempty(good_trials), warning('%s: no good trials, skipping.', subj); continue; end

    EEG = pop_loadset('filename', [key '_clean.set'], 'filepath', cfg.deriv_dir);
    ft_raw = eeglab2fieldtrip(EEG, 'raw', 'none');
    cfg_csd = []; cfg_csd.method = 'spline';
    ft_csd = ft_scalpcurrentdensity(cfg_csd, ft_raw);
    csd_data = ft_csd.trial{1};
    csd_labels = ft_csd.label;

    fr = local_chanidx(csd_labels, cfg.roi_frontal);
    pa = local_chanidx(csd_labels, cfg.roi_parietal);
    if numel(fr) < numel(cfg.roi_frontal) || numel(pa) < numel(cfg.roi_parietal)
        warning('%s: missing ROI channel(s) in CSD output, skipping.', subj); continue;
    end

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
                    X = zi.*conj(zj); im = imag(X);
                    s_im=sum(im,2); s_abs=sum(abs(im),2); s_sq=sum(im.^2,2);
                    dwpli_t = (s_im.^2 - s_sq) ./ (s_abs.^2 - s_sq);
                    plv_t = abs(mean(X./abs(X), 2));
                    DWPLI(su,cidx,p) = mean(dwpli_t,'omitnan');
                    PLV(su,cidx,p) = mean(plv_t,'omitnan');

                    Sxy = mean(mean(X)); Sxx = mean(mean(abs(zi).^2)); Syy = mean(mean(abs(zj).^2));
                    IMC(su,cidx,p) = imag(Sxy/sqrt(Sxx*Syy));
                end
            end
        end
    end
    fprintf('%s done.\n', subj);
end

% ======================= per-pair BF/TOST, all 3 metrics =======================
Delta = 0.4; r_prior = 0.707;
metrics = struct('dwPLI', DWPLI, 'PLV', PLV, 'ImC', IMC);
mnames = fieldnames(metrics);
rows = {};
d_val_ImC = [];   % kept separately for the cluster-mass step below

for mi = 1:numel(mnames)
    Mval = metrics.(mnames{mi});
    av = Mval(:,1,:); ai = Mval(:,2,:); sv = Mval(:,3,:); iv = Mval(:,4,:);
    d_val_pair = squeeze((av+sv)/2 - (ai+iv)/2);   % [subject x pair]
    if strcmp(mnames{mi}, 'ImC'), d_val_ImC = d_val_pair; end

    for p = 1:nF*nT
        d = d_val_pair(:,p); d = d(~isnan(d)); n = numel(d);
        if n < 3, continue; end
        dz = mean(d)/std(d); se = std(d)/sqrt(n); tobs = mean(d)/se; v = n-1;
        p_freq = 2*(1-tcdf(abs(tobs),v));
        BF10 = local_jzs_bf10(tobs, n, r_prior); BF01 = 1/BF10;
        bound_raw = Delta*std(d);
        t_lo=(mean(d)-(-bound_raw))/se; t_hi=(mean(d)-bound_raw)/se;
        p_tost = max(1-tcdf(t_lo,v), tcdf(t_hi,v));
        rows(end+1,:) = {labs{p}, mnames{mi}, n, dz, tobs, p_freq, BF10, BF01, p_tost, p_tost<0.05}; %#ok<SAGROW>
    end
end
Tpair = cell2table(rows, 'VariableNames', ...
    {'pair','metric','n','dz','t_freq','p_freq','BF10','BF01','p_TOST','equivalent_at_0p4'});
writetable(Tpair, fullfile(cfg.deriv_dir, 'csd_imc_perpair_valid_main_gamma_high.csv'));
fprintf('\nSaved csd_imc_perpair_valid_main_gamma_high.csv (%d rows)\n', height(Tpair));

Timc = Tpair(strcmp(Tpair.metric,'ImC'), :);
Timc = sortrows(Timc, 'p_freq');
fprintf('\n--- ImC per-pair, top 12 by uncorrected p_freq ---\n');
disp(Timc(1:min(12,height(Timc)), {'pair','n','dz','t_freq','p_freq','BF10','BF01'}));
n_nominal = sum(Timc.p_freq < 0.05);
fprintf('Pairs with p_freq<0.05 (uncorrected): %d of %d\n', n_nominal, height(Timc));

% ======================= cluster-mass test on ImC (real montage adjacency) =======================
fprintf('\n--- Cluster-mass permutation, ImC, valid_main/gamma_high, CSD-referenced ---\n');
neighbours_file = cfg.neighbours_file;
adj = local_build_montage_adjacency(cfg.roi_frontal, cfg.roi_parietal, neighbours_file);

D = d_val_ImC; bad = any(isnan(D),2); D = D(~bad,:); n = size(D,1);
tobs = mean(D) ./ (std(D)./sqrt(n));
alpha_form = 0.05; v = n-1; t_thr = tinv(1-alpha_form/2, v);
[obs_clusters, obs_mass] = local_find_clusters(tobs, adj, t_thr);
fprintf('n=%d, cluster-forming threshold |t|>%.3f, %d clusters found\n', n, t_thr, numel(obs_mass));

rng(1); nperm = cfg.stat.n_perm;
max_mass_null = zeros(nperm,1);
for pp = 1:nperm
    sgn = sign(rand(n,1)-0.5); Dp = D.*sgn; tp = mean(Dp)./(std(Dp)./sqrt(n));
    [~, mass_p] = local_find_clusters(tp, adj, t_thr);
    if isempty(mass_p), max_mass_null(pp)=0; else, max_mass_null(pp)=max(mass_p); end
end

clust_rows = {};
for k = 1:numel(obs_clusters)
    p_clust = mean(max_mass_null >= obs_mass(k));
    members = strjoin(labs(obs_clusters{k}), ', ');
    fprintf('  cluster %d: %d pairs, mass=%.2f, p_cluster=%.4f, members: %s\n', ...
        k, numel(obs_clusters{k}), obs_mass(k), p_clust, members);
    clust_rows(end+1,:) = {k, numel(obs_clusters{k}), obs_mass(k), p_clust, p_clust<0.05, members}; %#ok<SAGROW>
end
if isempty(obs_clusters)
    fprintf('  No supra-threshold clusters found.\n');
end
Tclust = cell2table(clust_rows, 'VariableNames', {'cluster_id','n_pairs','mass','p_cluster','sig','members'});
writetable(Tclust, fullfile(cfg.deriv_dir, 'csd_imc_clustermass_valid_main_gamma_high.csv'));
fprintf('\nSaved csd_imc_clustermass_valid_main_gamma_high.csv\n');

fprintf(['\nInterpretation guide:\n' ...
  '  - If p_freq<0.05 pairs (uncorrected) are scattered with no significant\n' ...
  '    cluster (p_cluster all >0.05): the ROI-average ImC anomaly is diffuse,\n' ...
  '    consistent with CSD noise amplification rather than a localized signal --\n' ...
  '    do not report it as a finding, note it only as a sensitivity-analysis\n' ...
  '    observation requiring replication.\n' ...
  '  - If a spatially coherent, significant cluster emerges (p_cluster<0.05):\n' ...
  '    this would be a genuinely new, CSD-revealed signal and warrants a\n' ...
  '    clearly-labeled exploratory subsection, NOT folded into the main\n' ...
  '    dissociation narrative without explicit caveats about its post-hoc,\n' ...
  '    single-dataset origin.\n']);

% ======================= local functions =======================
function idx = local_chanidx(labels, want)
    idx = [];
    for i = 1:numel(want)
        j = find(strcmpi(labels, want{i}), 1);
        if isempty(j), warning('channel "%s" not found', want{i}); else, idx(end+1) = j; end %#ok<AGROW>
    end
end

function bf10 = local_jzs_bf10(t, n, r)
    v = n-1; num = (1+t^2/v)^(-(v+1)/2);
    integrand = @(g) (1+n.*g).^(-0.5).*(1+t^2./((1+n.*g).*v)).^(-(v+1)/2).*(r/sqrt(2*pi)).*g.^(-1.5).*exp(-(r^2)./(2.*g));
    denom = integral(integrand,0,Inf); bf10 = 1/(num/denom);
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
