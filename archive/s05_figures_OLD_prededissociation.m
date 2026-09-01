% s05_figures.m  --  Publication-grade figures (300 dpi PNG), Pilot 2.
% Fig1: frontal_power by CUE (headline local effect, both bands significant)
% Fig2: parietal_power by CUE (gamma_low significant, gamma_high trend)
% Fig3: PLV by VALIDITY (gamma_high significant, gamma_low control/n.s.)
% Fig4: dwPLI null result -- by cue and by validity, both bands (headline
%       connectivity finding: local power/PLV effects present, but no
%       reliable phase-based frontoparietal connectivity effect)
% Each panel: per-subject paired lines, group mean +- SEM, permutation p.
% (Pair-level t-maps from s06_cluster_conn.m are separate files, already
% saved as fig3_conn_tmap_<contrast>.png -- not regenerated here.)

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
outdir = fullfile(cfg.deriv_dir,'figures'); if ~exist(outdir,'dir'), mkdir(outdir); end
subs = unique(R.subject, 'stable');
bands = fieldnames(cfg.bands);

getval = @(subj,bn,mt,cueLab,valLab) local_getval(R,subj,bn,mt,cueLab,valLab);

% ---- Figure 1: frontal_power by cue (both bands) ----
f1 = figure('Color','w','Position',[100 100 760 380]);
for bi = 1:numel(bands)
    [avatarV, stickV] = cue_and_valid_vectors(R, subs, bands{bi}, 'frontal_power');
    if cfg.log_power, avatarV = log(avatarV); stickV = log(stickV); ylab='log frontal power'; else, ylab='Frontal power'; end
    p = local_perm_p(avatarV - stickV, cfg.stat.n_perm, cfg.stat.tail);
    subplot(1,numel(bands),bi);
    local_pairplot(gca, stickV, avatarV, ylab, p, {'Stick','Avatar'});
    title(sprintf('%s band: frontal power by cue', bands{bi}), 'Interpreter','none');
end
exportgraphics(f1, fullfile(outdir,'fig1_frontal_power_by_cue.png'), 'Resolution', cfg.fig.dpi);

% ---- Figure 2: parietal_power by cue (both bands) ----
f2 = figure('Color','w','Position',[100 100 760 380]);
for bi = 1:numel(bands)
    [avatarV, stickV] = cue_and_valid_vectors(R, subs, bands{bi}, 'parietal_power');
    if cfg.log_power, avatarV = log(avatarV); stickV = log(stickV); ylab='log parietal power'; else, ylab='Parietal power'; end
    p = local_perm_p(avatarV - stickV, cfg.stat.n_perm, cfg.stat.tail);
    subplot(1,numel(bands),bi);
    local_pairplot(gca, stickV, avatarV, ylab, p, {'Stick','Avatar'});
    title(sprintf('%s band: parietal power by cue', bands{bi}), 'Interpreter','none');
end
exportgraphics(f2, fullfile(outdir,'fig2_parietal_power_by_cue.png'), 'Resolution', cfg.fig.dpi);

% ---- Figure 3: PLV by validity (both bands) ----
f3 = figure('Color','w','Position',[100 100 760 380]);
for bi = 1:numel(bands)
    [~,~,validV,invalidV] = cue_and_valid_vectors(R, subs, bands{bi}, 'plv');
    p = local_perm_p(validV - invalidV, cfg.stat.n_perm, cfg.stat.tail);
    subplot(1,numel(bands),bi);
    local_pairplot(gca, invalidV, validV, 'Frontoparietal PLV', p, {'Invalid','Valid'});
    title(sprintf('%s band: PLV by validity', bands{bi}), 'Interpreter','none');
end
exportgraphics(f3, fullfile(outdir,'fig3_plv_by_validity.png'), 'Resolution', cfg.fig.dpi);

% ---- Figure 4: dwPLI null result -- by cue and by validity, both bands ----
f4 = figure('Color','w','Position',[100 100 900 700]);
panel = 0;
for bi = 1:numel(bands)
    panel = panel+1;
    [avatarV, stickV] = cue_and_valid_vectors(R, subs, bands{bi}, 'dwpli');
    p = local_perm_p(avatarV - stickV, cfg.stat.n_perm, cfg.stat.tail);
    subplot(2,2,panel);
    local_pairplot(gca, stickV, avatarV, 'Frontoparietal dwPLI', p, {'Stick','Avatar'});
    title(sprintf('%s band: dwPLI by cue (n.s.)', bands{bi}), 'Interpreter','none');
end
for bi = 1:numel(bands)
    panel = panel+1;
    [~,~,validV,invalidV] = cue_and_valid_vectors(R, subs, bands{bi}, 'dwpli');
    p = local_perm_p(validV - invalidV, cfg.stat.n_perm, cfg.stat.tail);
    subplot(2,2,panel);
    local_pairplot(gca, invalidV, validV, 'Frontoparietal dwPLI', p, {'Invalid','Valid'});
    title(sprintf('%s band: dwPLI by validity (n.s.)', bands{bi}), 'Interpreter','none');
end
exportgraphics(f4, fullfile(outdir,'fig4_dwpli_null.png'), 'Resolution', cfg.fig.dpi);

fprintf('Figures saved to %s\n', outdir);

% ======================= local functions =======================
function [avatarV, stickV, validV, invalidV] = cue_and_valid_vectors(R, subs, bn, mt)
    % per-subject cue-averaged (collapsed over validity) and
    % validity-averaged (collapsed over cue) values for a given band/metric
    n = numel(subs);
    avatarV = nan(n,1); stickV = nan(n,1); validV = nan(n,1); invalidV = nan(n,1);
    for i = 1:n
        av = local_getval(R, subs{i}, bn, mt, 'avatar','valid');
        ai = local_getval(R, subs{i}, bn, mt, 'avatar','invalid');
        sv = local_getval(R, subs{i}, bn, mt, 'stick','valid');
        si = local_getval(R, subs{i}, bn, mt, 'stick','invalid');
        avatarV(i) = mean([av ai],'omitnan');
        stickV(i)  = mean([sv si],'omitnan');
        validV(i)  = mean([av sv],'omitnan');
        invalidV(i)= mean([ai si],'omitnan');
    end
end

function v = local_getval(R, subj, bn, mt, cueLab, valLab)
    r = R(strcmp(R.subject,subj) & strcmp(R.band,bn) & ...
          strcmp(R.cue,cueLab) & strcmp(R.valid,valLab), :);
    if isempty(r), v = NaN; else, v = r.(mt)(1); end
end

function p = local_perm_p(d, nperm, tail)
    d = d(~isnan(d)); n = numel(d);
    tobs = mean(d)/(std(d)/sqrt(n));
    tnull = zeros(nperm,1);
    for pp = 1:nperm
        dp = d .* sign(rand(n,1)-0.5);
        tnull(pp) = mean(dp)/(std(dp)/sqrt(n));
    end
    switch tail
      case 0,  p = mean(abs(tnull) >= abs(tobs));
      case 1,  p = mean(tnull >= tobs);
      case -1, p = mean(tnull <= tobs);
    end
end

function local_pairplot(ax, lo, hi, ylab, pval, xl)
    lo = lo(:); hi = hi(:); ok = ~isnan(lo) & ~isnan(hi);
    lo = lo(ok); hi = hi(ok); n = numel(lo);
    hold(ax,'on');
    for i = 1:n
        plot(ax, [1 2], [lo(i) hi(i)], '-', 'Color', [.7 .7 .7]);
    end
    plot(ax, ones(n,1), lo, 'o', 'MarkerFaceColor',[.3 .5 .8], 'MarkerEdgeColor','none');
    plot(ax, 2*ones(n,1), hi, 'o', 'MarkerFaceColor',[.85 .4 .3], 'MarkerEdgeColor','none');
    m = [mean(lo) mean(hi)]; e = [std(lo) std(hi)]/sqrt(n);
    errorbar(ax, [1 2], m, e, 'k', 'LineWidth', 2, 'CapSize', 12);
    xlim(ax,[.6 2.4]); set(ax,'XTick',[1 2],'XTickLabel',xl);
    ylabel(ax, ylab); box(ax,'off');
    yl = ylim(ax); text(ax, 1.5, yl(2), sprintf('p = %.3f', pval), ...
        'HorizontalAlignment','center','VerticalAlignment','top');
end
