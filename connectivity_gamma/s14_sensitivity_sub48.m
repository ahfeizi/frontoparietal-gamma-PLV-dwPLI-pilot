% s14_sensitivity_sub48.m -- Sensitivity check: does excluding sub-48
% (manually flagged bad channel P1, cfg.manual_bad_channels in config.m)
% change any of the established results?
%
% Cheap by design: reuses conn_results.mat (already has all subjects'
% dwpli_pairs/plv_pairs) and trials.mat -- no raw signal reconstruction
% needed here, unlike s15 (RT_stat sensitivity).
%
% Checks, with vs without sub-48:
%   (A) The gamma Bayes/TOST table (3 contrasts x 2 metrics x 2 subbands),
%       same method as s10_bayes_tost_gamma.m.
%   (B) The valid_main / PLV / gamma_high cluster-mass finding (F2-P4 and
%       surrounding cluster), same method as s06c_robustness (real montage
%       adjacency from neighbours.mat).
%   (C) The Phase 0 behavioral ANOVA on acc_loc (cue x validity
%       interaction), same method as s00b_phase0_behavior.m.

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
M = load(fullfile(cfg.deriv_dir,'manifest.mat')); manifest = M.manifest;

exclude_subj = 'sub-48';
Delta = 0.4; r_prior = 0.707;   % copied from s10/s12/s13, unchanged

fprintf('=== Sensitivity: excluding %s ===\n\n', exclude_subj);

% ======================= (A) gamma Bayes/TOST table, with vs without sub-48 =======================
metrics = {'dwpli','plv'};
subbands = {'gamma_low','gamma_high'};
getval = @(Ruse,subj,bn,mt,cueLab,valLab) local_getval(Ruse, subj, bn, mt, cueLab, valLab);

rows_A = {};   % <-- NEW

for excl = [false true]
    if excl
        Ruse = R(~strcmp(R.subject, exclude_subj), :);
        tag = sprintf('WITHOUT %s', exclude_subj);
    else
        Ruse = R;
        tag = 'FULL (all 47)';
    end
    subs = unique(Ruse.subject);
    fprintf('--- (A) Gamma Bayes/TOST, %s (n=%d) ---\n', tag, numel(subs));
    fprintf('%-8s %-8s %-12s %3s  %7s  %7s  %8s  %8s  %8s\n', 'contrast','metric','band','n','dz','t','BF10','BF01','p_TOST');
    for s = 1:numel(subbands)
        bn = subbands{s};
        for m = 1:numel(metrics)
            mt = metrics{m};
            av=nan(numel(subs),1); ai=nan(numel(subs),1); sv=nan(numel(subs),1); si_=nan(numel(subs),1);
            for i=1:numel(subs)
                av(i)=getval(Ruse,subs{i},bn,mt,'avatar','valid');
                ai(i)=getval(Ruse,subs{i},bn,mt,'avatar','invalid');
                sv(i)=getval(Ruse,subs{i},bn,mt,'stick','valid');
                si_(i)=getval(Ruse,subs{i},bn,mt,'stick','invalid');
            end
            contrasts = struct('d_cue',(av+ai)/2-(sv+si_)/2, 'd_val',(av+sv)/2-(ai+si_)/2, 'd_int',(av-ai)-(sv-si_));
            cn = fieldnames(contrasts);
            for c = 1:numel(cn)
                d = contrasts.(cn{c}); d = d(~isnan(d));
                r0 = local_eval(d, Delta, r_prior);
                fprintf('%-8s %-8s %-12s %3d  %7.3f  %7.2f  %8.3f  %8.3f  %8.4f\n', ...
                    cn{c}, mt, bn, r0{1}, r0{2}, r0{3}, r0{5}, r0{6}, r0{7});
                rows_A(end+1,:) = {tag, cn{c}, mt, bn, r0{1}, r0{2}, r0{3}, r0{4}, r0{5}, r0{6}, r0{7}, r0{8}}; %#ok<SAGROW>  % <-- NEW
            end
        end
    end
    fprintf('\n');
end

% <-- NEW
Tout_A = cell2table(rows_A, 'VariableNames', ...
    {'policy','contrast','metric','band','n','dz','t','p_freq','BF10','BF01','p_TOST','equivalent_at_0p4'});
writetable(Tout_A, fullfile(cfg.deriv_dir, 'sensitivity_gammahigh_bayes.csv'));
fprintf('Saved sensitivity_gammahigh_bayes.csv\n\n');

% ======================= (B) F2-P4 / cluster-mass, with vs without sub-48 =======================
fprintf('--- (B) valid_main / PLV / gamma_high cluster-mass, with vs without %s ---\n', exclude_subj);
nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_parietal);
labs = cell(1,nF*nT);
for i=1:nF, for j=1:nT, labs{(i-1)*nT+j} = sprintf('%s-%s', cfg.roi_frontal{i}, cfg.roi_parietal{j}); end, end
target_pair = find(strcmp(labs,'F2-P4'));
neighbours_file = cfg.neighbours_file;
adj = local_build_montage_adjacency(cfg.roi_frontal, cfg.roi_parietal, neighbours_file);

rows_B = {};   % <-- NEW

for excl = [false true]
    if excl
        Ruse = R(~strcmp(R.subject, exclude_subj), :);
        tag = sprintf('WITHOUT %s', exclude_subj);
    else
        Ruse = R; tag = 'FULL (all 47)';
    end
    subs = unique(Ruse.subject, 'stable');
    bn = 'gamma_high';
    D = nan(numel(subs), nF*nT);
    for si = 1:numel(subs)
        av = local_getpairs(Ruse, subs{si}, bn, 'avatar','valid');
        ai = local_getpairs(Ruse, subs{si}, bn, 'avatar','invalid');
        sv = local_getpairs(Ruse, subs{si}, bn, 'stick','valid');
        iv = local_getpairs(Ruse, subs{si}, bn, 'stick','invalid');
        if isempty(av)||isempty(ai)||isempty(sv)||isempty(iv), continue; end
        D(si,:) = (av+sv)/2 - (ai+iv)/2;
    end
    bad = any(isnan(D),2); D = D(~bad,:); n = size(D,1);
    tobs = mean(D) ./ (std(D)./sqrt(n));
    fprintf('%s (n=%d): F2-P4 t=%.3f\n', tag, n, tobs(target_pair));

    rng(1); nperm = cfg.stat.n_perm;
    alpha_form = 0.05; v = n-1; t_thr = tinv(1-alpha_form/2, v);
    [obs_clusters, obs_mass] = local_find_clusters(tobs, adj, t_thr);
    max_mass_null = zeros(nperm,1);
    for p = 1:nperm
        sgn = sign(rand(n,1)-0.5); Dp = D.*sgn; tp = mean(Dp)./(std(Dp)./sqrt(n));
        [~, mass_p] = local_find_clusters(tp, adj, t_thr);
        if isempty(mass_p), max_mass_null(p)=0; else, max_mass_null(p)=max(mass_p); end
    end
    % <-- NEW: defaults, in case F2-P4 is not part of any cluster this run
    f2p4_t = tobs(target_pair); n_pairs_out = NaN; mass_out = NaN; p_cluster_out = NaN; in_cluster = false;
    if any(cellfun(@(c) any(c==target_pair), obs_clusters))
        idx = find(cellfun(@(c) any(c==target_pair), obs_clusters));
        m = obs_mass(idx); p_cluster = mean(max_mass_null >= m);
        fprintf('  F2-P4''s cluster: %d pairs, mass=%.2f, p_cluster=%.4f\n', numel(obs_clusters{idx}), m, p_cluster);
        n_pairs_out = numel(obs_clusters{idx}); mass_out = m; p_cluster_out = p_cluster; in_cluster = true;  % <-- NEW
    else
        fprintf('  F2-P4 not part of any supra-threshold cluster in this run.\n');
    end
    rows_B(end+1,:) = {tag, n, f2p4_t, in_cluster, n_pairs_out, mass_out, p_cluster_out}; %#ok<SAGROW>  % <-- NEW
end
fprintf('\n');

% <-- NEW
Tout_B = cell2table(rows_B, 'VariableNames', ...
    {'policy','n','f2p4_t','f2p4_in_cluster','cluster_n_pairs','cluster_mass','p_cluster'});
writetable(Tout_B, fullfile(cfg.deriv_dir, 'sensitivity_f2p4_cluster.csv'));
fprintf('Saved sensitivity_f2p4_cluster.csv\n\n');

% ======================= (C) Phase 0 acc_loc ANOVA, with vs without sub-48 =======================
fprintf('--- (C) acc_loc ANOVA (cue x validity interaction), with vs without %s ---\n', exclude_subj);
cue_code = struct('avatar', cfg.ti.cue_avatar, 'stick', cfg.ti.cue_stick);
valid_code = struct('valid', cfg.ti.val_valid, 'invalid', cfg.ti.val_invalid);
cueLabs = {'avatar','stick'}; valLabs = {'valid','invalid'};

acc_rows = {};
for su = 1:height(manifest)
    subj = char(manifest.subject(su)); key = char(manifest.key(su));
    trials_file = fullfile(cfg.deriv_dir, [key '_trials.mat']);
    if ~isfile(trials_file), continue; end
    Tt = load(trials_file); trials = Tt.trials;
    cueArr=[trials.cue]; validArr=[trials.valid]; accLocArr=[trials.acc_loc];
    for c=1:2, for v=1:2
        sel = (cueArr==cue_code.(cueLabs{c})) & (validArr==valid_code.(valLabs{v}));
        acc_rows(end+1,:) = {subj, cueLabs{c}, valLabs{v}, mean(accLocArr(sel),'omitnan')}; %#ok<SAGROW>
    end, end
end
Tal = cell2table(acc_rows, 'VariableNames', {'subject','cue','valid','acc'});

for excl = [false true]
    if excl
        Tuse = Tal(~strcmp(Tal.subject, exclude_subj), :);
        tag = sprintf('WITHOUT %s', exclude_subj);
    else
        Tuse = Tal; tag = 'FULL (all 47)';
    end
    fprintf('%s (n=%d subjects):\n', tag, numel(unique(Tuse.subject)));
    local_run_anova(Tuse, sprintf('acc_loc_%s', strrep(tag,' ','_')), cfg.deriv_dir);
end

% ======================= local functions =======================
function v = local_getval(R, subj, bn, mt, cueLab, valLab)
    r = R(strcmp(R.subject,subj) & strcmp(R.band,bn) & strcmp(R.cue,cueLab) & strcmp(R.valid,valLab), :);
    if isempty(r), v = NaN; else, v = r.(mt)(1); end
end

function pairs = local_getpairs(R, subj, bn, cueLab, valLab)
    r = R(strcmp(R.subject,subj) & strcmp(R.band,bn) & strcmp(R.cue,cueLab) & strcmp(R.valid,valLab), :);
    if isempty(r), pairs = []; else, pairs = r.plv_pairs{1}; if iscell(pairs), pairs=pairs{1}; end; end
end

function row = local_eval(d, Delta, r_prior)
    n = numel(d); dz = mean(d)/std(d); se = std(d)/sqrt(n); tobs = mean(d)/se; v = n-1;
    p_freq = 2*(1-tcdf(abs(tobs),v));
    BF10 = local_jzs_bf10(tobs, n, r_prior); BF01 = 1/BF10;
    bound_raw = Delta*std(d);
    t_lo=(mean(d)-(-bound_raw))/se; t_hi=(mean(d)-bound_raw)/se;
    p_tost = max(1-tcdf(t_lo,v), tcdf(t_hi,v));
    row = {n, dz, tobs, p_freq, BF10, BF01, p_tost, p_tost<0.05};
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

function local_run_anova(Tlong, label, deriv_dir)
    subs = unique(Tlong.subject,'stable');
    wide = table(subs,'VariableNames',{'subject'});
    condnames={'avatar_valid','avatar_invalid','stick_valid','stick_invalid'};
    combos={'avatar','valid';'avatar','invalid';'stick','valid';'stick','invalid'};
    for k=1:4
        v=nan(numel(subs),1);
        for i=1:numel(subs)
            r=Tlong(strcmp(Tlong.subject,subs{i})&strcmp(Tlong.cue,combos{k,1})&strcmp(Tlong.valid,combos{k,2}),:);
            if ~isempty(r), v(i)=r.acc(1); end
        end
        wide.(condnames{k})=v;
    end
    bad=any(isnan(wide{:,2:5}),2); wide=wide(~bad,:);
    within=table(categorical({'avatar';'avatar';'stick';'stick'}), categorical({'valid';'invalid';'valid';'invalid'}), 'VariableNames',{'cue','valid'});
    rm=fitrm(wide, sprintf('%s-%s ~ 1',condnames{1},condnames{4}), 'WithinDesign', within);
    ranovatbl=ranova(rm,'WithinModel','cue*valid');
    intRow = ranovatbl('(Intercept):cue:valid',:);
    F_int = double(intRow.F); p_int = double(intRow.pValue);
    fprintf('  interaction: F=%.3f, p=%.4f\n', F_int, p_int);

    plainT = table(string(ranovatbl.Properties.RowNames), 'VariableNames', {'term'});
    varnames = ranovatbl.Properties.VariableNames;
    for vn = 1:numel(varnames)
        col = ranovatbl.(varnames{vn});
        try, col = double(col); catch, end
        plainT.(varnames{vn}) = col;
    end
    writetable(plainT, fullfile(deriv_dir, sprintf('sensitivity_%s.csv', label)));
end
