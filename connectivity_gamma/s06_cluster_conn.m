% s06_cluster_conn.m  --  FWER-controlled test across the full connectivity
% space (frontal x parietal pairs x bands), Pilot 2, for BOTH dwPLI and PLV.
%
% UPDATED from the original dwPLI-only version: this script now loops over
% metric in {'dwpli','plv'} as well as contrast and band, using the same
% max-statistic sign-flip permutation logic (Maris-Oostenveld family),
% unchanged. It was extended after S10 (s10_bayes_tost_gamma.m) found a
% strong ROI-averaged Bayesian effect for d_val/PLV/gamma_high (BF10=45)
% that had no counterpart in dwPLI -- the original script only read
% r.dwpli_pairs{1} and could not test whether that PLV effect was localized
% to specific frontal-parietal pairs or diffuse across the pair space.
%
% BACKWARD COMPATIBILITY: the original per-contrast dwPLI-only file
%   conn_cluster_results_<contrast>.csv
% is still written, with IDENTICAL content/columns to before (no metric
% column), since other scripts/reports (e.g. the Pilot 2 write-up) already
% reference this filename for dwPLI. A new combined file
%   conn_cluster_results_<contrast>_by_metric.csv
% is written alongside it, stacking both metrics (256 rows: 128 dwPLI +
% 128 PLV) with an explicit 'metric' column. Figures are also split by
% metric: fig3_conn_tmap_<contrast>_dwpli.png / _plv.png (the original
% fig3_conn_tmap_<contrast>.png is kept for dwPLI to preserve compatibility).

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
rng(1);

bands = fieldnames(cfg.bands);
nF = numel(cfg.roi_frontal); nT = numel(cfg.roi_parietal); npair = nF*nT;
subs = unique(R.subject, 'stable');

metrics = {'dwpli','plv'};
pair_field = struct('dwpli','dwpli_pairs', 'plv','plv_pairs');

getpairs = @(subj,bn,cueLab,valLab,fieldname) local_getpairs(R, subj, bn, cueLab, valLab, fieldname);

contrast_names = {'cue_main','valid_main','interaction'};

for ci = 1:numel(contrast_names)
    cname = contrast_names{ci};
    fprintf('\n========== Contrast: %s ==========\n', cname);

    combined_rows = {};  % accumulates both metrics for the _by_metric.csv

    for mi = 1:numel(metrics)
        mt = metrics{mi};
        fld = pair_field.(mt);
        fprintf('\n---- metric: %s ----\n', mt);

        labs = {}; D = [];
        for bi = 1:numel(bands)
            bn = bands{bi};
            Db = nan(numel(subs), npair);
            for si = 1:numel(subs)
                subj = subs{si};
                av = getpairs(subj, bn, 'avatar', 'valid',   fld);
                ai = getpairs(subj, bn, 'avatar', 'invalid', fld);
                sv = getpairs(subj, bn, 'stick',  'valid',   fld);
                si_ = getpairs(subj, bn, 'stick',  'invalid', fld);
                if isempty(av) || isempty(ai) || isempty(sv) || isempty(si_)
                    continue;   % leaves this subject's row as NaN -> dropped below
                end
                switch cname
                    case 'cue_main'
                        d = (av+ai)/2 - (sv+si_)/2;
                    case 'valid_main'
                        d = (av+sv)/2 - (ai+si_)/2;
                    case 'interaction'
                        d = (av-ai) - (sv-si_);
                end
                Db(si,:) = d;
            end
            D = [D Db]; %#ok<AGROW>
            for i = 1:nF, for j = 1:nT
                labs{end+1} = sprintf('%-9s %s-%s', bn, cfg.roi_frontal{i}, cfg.roi_parietal{j}); %#ok<AGROW,SAGROW>
            end, end %#ok<ALIGN>
        end

        bad = any(isnan(D),2); nDrop = sum(bad); D = D(~bad,:); n = size(D,1); Mcol = size(D,2);
        fprintf('Subjects used: %d (dropped %d with incomplete data)\n', n, nDrop);
        if n < 3
            fprintf('Too few subjects for %s / %s, skipping.\n', cname, mt);
            continue;
        end

        % ----- observed and permutation statistics (max-stat, sign-flip) -----
        tobs = mean(D) ./ (std(D)./sqrt(n));
        nperm = cfg.stat.n_perm;
        tnull = zeros(nperm, Mcol); maxnull = zeros(nperm,1);
        for p = 1:nperm
            sgn = sign(rand(n,1)-0.5);
            Dp = D .* sgn;                          % same flip across all columns
            t = mean(Dp) ./ (std(Dp)./sqrt(n));
            tnull(p,:) = t; maxnull(p) = max(abs(t));
        end
        thr     = prctile(maxnull, 95);
        p_unc   = mean(abs(tnull) >= abs(tobs), 1);
        p_fwer  = mean(maxnull   >= abs(tobs), 1);
        sig     = p_fwer < 0.05;

        T = table(labs(:), tobs(:), p_unc(:), p_fwer(:), sig(:), ...
            'VariableNames', {'pair_band','t','p_uncorr','p_fwer','sig_fwer'});
        T = sortrows(T, 't', 'descend', 'ComparisonMethod','abs');
        fprintf('FWER threshold |t| (95%% of max-null) = %.2f\n', thr);
        fprintf('FWER-significant pairs: %d of %d\n\n', sum(sig), Mcol);
        disp(T(1:min(12,height(T)), :));

        if strcmp(mt,'dwpli')
            % unchanged filename/content for backward compatibility
            writetable(T, fullfile(cfg.deriv_dir, sprintf('conn_cluster_results_%s.csv', cname)));
        end

        % accumulate for combined by-metric output (unsorted, original pair order)
        Traw = table(labs(:), repmat({mt},numel(labs),1), tobs(:), p_unc(:), p_fwer(:), sig(:), ...
            'VariableNames', {'pair_band','metric','t','p_uncorr','p_fwer','sig_fwer'});
        combined_rows = [combined_rows; table2cell(Traw)]; %#ok<AGROW>

        % ----- t-map figure (frontal x parietal), one panel per band -----
        f = figure('Color','w','Position',[100 100 820 360]);
        for bi = 1:numel(bands)
            cols = (bi-1)*npair + (1:npair);
            Tm = reshape(tobs(cols), nT, nF)';
            Sm = reshape(sig(cols),  nT, nF)';
            subplot(1,numel(bands),bi);
            imagesc(Tm); axis image; colorbar; clim([-max(abs(tobs)) max(abs(tobs))]);
            colormap(parula);
            set(gca,'XTick',1:nT,'XTickLabel',cfg.roi_parietal, ...
                    'YTick',1:nF,'YTickLabel',cfg.roi_frontal);
            title(sprintf('%s band: %s (%s, t)', bands{bi}, cname, mt), 'Interpreter','none');
            [ri,cci] = find(Sm);
            hold on; plot(cci, ri, 'k*', 'MarkerSize',10, 'LineWidth',1.5);
        end
        outdir = fullfile(cfg.deriv_dir,'figures'); if ~exist(outdir,'dir'), mkdir(outdir); end
        if strcmp(mt,'dwpli')
            % unchanged filename for backward compatibility
            exportgraphics(f, fullfile(outdir, sprintf('fig3_conn_tmap_%s.png', cname)), 'Resolution', cfg.fig.dpi);
        end
        exportgraphics(f, fullfile(outdir, sprintf('fig3_conn_tmap_%s_%s.png', cname, mt)), 'Resolution', cfg.fig.dpi);
        fprintf('Saved figures for %s / %s\n', cname, mt);
    end

    Tcombined = cell2table(combined_rows, ...
        'VariableNames', {'pair_band','metric','t','p_uncorr','p_fwer','sig_fwer'});
    writetable(Tcombined, fullfile(cfg.deriv_dir, sprintf('conn_cluster_results_%s_by_metric.csv', cname)));
    fprintf('\nSaved conn_cluster_results_%s.csv (dwPLI, unchanged) and conn_cluster_results_%s_by_metric.csv (dwPLI+PLV, 256 rows)\n', cname, cname);
end

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
