% s05_figures.m  --  Publication-grade figures (300 dpi PNG), Pilot 2.
% REWRITTEN: the original version of this script (now in /archive/) plotted
% PLV and dwPLI as separate, independent paired-comparison panels. That
% framing predates the PLV-dwPLI dissociation finding and does not convey
% the paper's actual claim, which rests on MULTI-METRIC CONVERGENCE, not any
% single measure in isolation (see README "Honest limitations"). This
% version reads directly from the verified files in /results/ (not from
% conn_results.mat) so the figures are guaranteed to match the numbers
% already cross-checked against the preprint.
%
% Fig3: Bayes factors (BF10/BF01) across PLV, dwPLI, and imaginary coherency
%       for the validity x high-gamma contrast -- the paper's central figure.
% Fig4: Circular distribution of F2-P4 phase differences (47 subjects x 4
%       conditions = 188 cells), the trial-level evidence for a zero-lag
%       (volume-conduction) signature.
% Fig5: Stability of the PLV BF10 across the primary analysis and two
%       sensitivity checks (sub-48 removed; RT_stat reactivated).
%
% Fig1 (behavioral replication) and Fig2 (cluster schematic) are NOT
% generated here: Fig1 needs the original study's F-statistics as reference
% points (present in the preprint text but not in any results/ file yet --
% add them to bayes_tost_gamma.csv's sibling behavioral table, or hardcode
% with a clear citation, before writing that script). Fig2 is a schematic
% network diagram, not a data plot, and is better done by hand or in the
% Visualizer than regenerated from cluster_mass_*_members.csv coordinates.

clear; clc;
cfg = config();
resdir = cfg.deriv_dir;   % assumes /results/ files live alongside other derivatives;
                           % EDIT if you move results/ to a separate path.
outdir = fullfile(cfg.deriv_dir,'figures'); if ~exist(outdir,'dir'), mkdir(outdir); end

% ============================================================
% Figure 3: Bayes factors across PLV / dwPLI / imaginary coherency
% ============================================================
BT = readtable(fullfile(resdir,'bayes_tost_gamma.csv'));
plv_row   = BT(strcmp(BT.contrast,'d_val') & strcmp(BT.metric,'plv')   & strcmp(BT.band,'gamma_high'), :);
dwpli_row = BT(strcmp(BT.contrast,'d_val') & strcmp(BT.metric,'dwpli') & strcmp(BT.band,'gamma_high'), :);
BF10_plv   = plv_row.BF10(1);      % ROI-averaged, ~45.16
BF01_dwpli = dwpli_row.BF01(1);    % ROI-averaged, ~6.13

% Imaginary coherency: preprint reports BOTH an ROI-averaged BF01 (~6.18)
% and the F2-P4-specific BF01 (~6.31). Only the F2-P4 value is currently in
% results/ (imag_coherence_valid_main_gamma_high.csv is per-pair). Using
% F2-P4 here as the closest verified number -- REPLACE with the ROI-average
% once that summary file exists (should look like bayes_tost_gamma.csv's
% ROI-level rows; ask for it before treating this bar as final for
% publication, since the printed BF value will read 6.31, not the 6.18 the
% preprint text states for the ROI-averaged bar).
IC = readtable(fullfile(resdir,'imag_coherence_valid_main_gamma_high.csv'));
ic_row = IC(strcmp(IC.pair,'F2-P4'), :);
BF01_ImC_F2P4 = ic_row.BF01(1);    % ~6.31 (F2-P4 specific, NOT ROI-average)

labels = {'PLV','dwPLI','ImC (F2-P4)'};
% Plot on a common "evidence for effect" axis: BF10 as-is, BF01 metrics
% flipped to 1/BF01 so all three bars share the same "evidence for H1"
% direction, matching the log-scale framing of the preprint's Figure 3.
vals = [BF10_plv, 1/BF01_dwpli, 1/BF01_ImC_F2P4];

f3 = figure('Color','w','Position',[100 100 520 420]);
b = bar(vals, 'FaceColor','flat');
b.CData = [0.85 0.33 0.10; 0.20 0.55 0.75; 0.20 0.55 0.75];  % PLV highlighted, robust metrics matched
set(gca,'YScale','log','XTick',1:3,'XTickLabel',labels);
ylabel('Evidence for effect (BF10 scale; dwPLI/ImC shown as 1/BF01)');
yline(3, 'k--', 'BF = 3 ("substantial evidence")', 'LabelHorizontalAlignment','left');
yline(1/3, 'k:', 'BF01 = 3', 'LabelHorizontalAlignment','left');
title({'Validity \times high-gamma: PLV vs. volume-conduction-robust metrics', ...
       sprintf('PLV BF10=%.2f | dwPLI BF01=%.2f | ImC(F2-P4) BF01=%.2f', BF10_plv, BF01_dwpli, BF01_ImC_F2P4)});
box off;
exportgraphics(f3, fullfile(outdir,'fig3_bayes_factors_convergence.png'), 'Resolution', cfg.fig.dpi);

% ============================================================
% Figure 4: circular distribution of F2-P4 phase differences
% ============================================================
PD = readtable(fullfile(resdir,'phase_diff_diagnostic_F2P4_gamma_high.csv'));
theta = deg2rad(PD.mean_angle_deg);
rho   = PD.resultant_len;

f4 = figure('Color','w','Position',[100 100 520 520]);
pax = polaraxes; hold(pax,'on');
polarscatter(pax, theta, rho, 24, 'filled', 'MarkerFaceAlpha', 0.55, ...
    'MarkerFaceColor',[0.20 0.35 0.65]);
% shaded band: within 20 deg of 180
th_band = deg2rad(linspace(160,200,100));
polarplot(pax, th_band, max(rho)*ones(size(th_band)), 'r-', 'LineWidth', 4);
pax.ThetaZeroLocation = 'right';
n_in_band = sum(abs(180 - abs(PD.mean_angle_deg)) <= 20);
title(pax, sprintf('F2-P4 phase differences, N=%d cells\n%d/%d (%.1f%%) within 20\\circ of 180\\circ', ...
    height(PD), n_in_band, height(PD), 100*n_in_band/height(PD)));
exportgraphics(f4, fullfile(outdir,'fig4_phase_diff_circular.png'), 'Resolution', cfg.fig.dpi);

% ============================================================
% Figure 5: stability of PLV BF10 across primary + sensitivity analyses
% ============================================================
SB = readtable(fullfile(resdir,'sensitivity_gammahigh_bayes.csv'));
SR = readtable(fullfile(resdir,'sensitivity_rtstat_gammahigh.csv'));

bf_primary = BF10_plv;  % from bayes_tost_gamma.csv, = FULL (all 47) in SB too
bf_sub48   = SB.BF10(strcmp(SB.policy,'WITHOUT sub-48') & strcmp(SB.contrast,'d_val') & strcmp(SB.metric,'plv'));
bf_rtstat  = SR.BF10(strcmp(SR.policy,'RT_stat filter DISABLED for sub-01..04') & strcmp(SR.contrast,'d_val') & strcmp(SR.metric,'plv'));

sens_labels = {'Primary (N=47)','sub-48 removed','RT\_stat reactivated'};
sens_vals = [bf_primary, bf_sub48, bf_rtstat];

f5 = figure('Color','w','Position',[100 100 520 380]);
bar(sens_vals, 'FaceColor',[0.85 0.33 0.10]);
set(gca,'XTick',1:3,'XTickLabel',sens_labels);
ylabel('PLV BF10 (validity \times high-gamma)');
ylim([0, max(sens_vals)*1.25]);
for i = 1:3
    text(i, sens_vals(i), sprintf('%.2f', sens_vals(i)), ...
        'HorizontalAlignment','center','VerticalAlignment','bottom');
end
title('PLV effect is stable across sensitivity checks');
box off;
exportgraphics(f5, fullfile(outdir,'fig5_sensitivity_stability.png'), 'Resolution', cfg.fig.dpi);

fprintf('Saved fig3_bayes_factors_convergence.png, fig4_phase_diff_circular.png, fig5_sensitivity_stability.png to %s\n', outdir);
fprintf('NOTE: Fig1 (behavioral) and Fig2 (cluster schematic) are not generated by this script -- see header comment.\n');
fprintf('NOTE: Fig3''s ImC bar uses the F2-P4 value (BF01=%.2f), not the ROI-averaged value (~6.18) the preprint text cites -- see inline comment.\n', BF01_ImC_F2P4);
