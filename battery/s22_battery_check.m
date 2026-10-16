% s22_battery_check.m -- Regression check of s21_battery.m output against
% the existing pipeline results. Reads battery_persubject.csv and
% battery_phasediff.csv from cfg.battery.check_dir and prints PASS/FAIL per check.
%
%   T1  Hilbert / average-ref ROI-average PLV and dwPLI per subject x band x
%       cue x validity == conn_results.mat (s03)                 tol 1e-6
%   T2  Hilbert / average-ref ImC, gamma_high, validity contrast
%       d_val = (avatar_valid+stick_valid)/2 - (avatar_invalid+stick_invalid)/2,
%       ROI average over pairs: dz == imag_coherence_roi_average_gamma_high.csv
%       (s13)                                                     tol 0.01
%   T3  Hilbert / CSD, gamma_high, same d_val contrast, dwPLI/PLV/ImC at
%       roi_avg and F2-P4: dz == csd_sensitivity_summary.csv (s18)  tol 0.01
%   T4  Hilbert / average-ref, gamma_high, F2-P4 phase difference: 175 of
%       188 subject x condition cells within +-20 deg of 180 deg (s13)
%   C1  battery_phasediff.csv has exactly one row per subject x condition x
%       band x decomposition x reference x pair of battery_persubject.csv,
%       with resultant_len in [0, 1]
%
% All inputs are read from cfg.battery.check_dir (default: <repo>/results/,
% the released copies; set it to cfg.deriv_dir to check a fresh s21 run).
% For each reference file the full path and last-modified date are printed;
% if check_dir is not <repo>/results/ and a file of the same name exists
% there, the two are compared and FAIL is printed if their contents differ
% (CSV line endings normalised).

clear; clc;
repo = fileparts(fileparts(mfilename('fullpath')));   % repository root (this script is in battery/)
addpath(repo, fullfile(repo, 'battery'));             % config.m, util_bandpass.m, util_chanidx.m; util_conn_pairs.m
cfg = config();
tol_val = 1e-6;
tol_dz  = 0.01;
phase_expect_within = 175;   % cells within +-20 deg of 180 deg
phase_expect_cells  = 188;   % 47 subjects x 4 conditions
phase_tol_deg       = 20;
nfail = 0;

conds = {'avatar_valid','avatar_invalid','stick_valid','stick_invalid'};
nF = numel(cfg.roi_frontal); nP = numel(cfg.roi_parietal);
pairs = cell(1, nF*nP);
for i = 1:nF, for j = 1:nP, pairs{(i-1)*nP+j} = sprintf('%s-%s', cfg.roi_frontal{i}, cfg.roi_parietal{j}); end, end
target_pair = find(strcmp(pairs, 'F2-P4'));

fprintf('Reading battery output from %s ...\n', cfg.battery.check_dir);
B = local_read(fullfile(cfg.battery.check_dir,'battery_persubject.csv'), ...
    {'subject','condition','band','measure','decomposition','reference','pair'}, {'value'});
P = local_read(fullfile(cfg.battery.check_dir,'battery_phasediff.csv'), ...
    {'subject','condition','band','decomposition','reference','pair'}, {'mean_angle_deg','resultant_len'});

% ---- T1 ----
[L, nfail] = local_ref(cfg, repo, 'conn_results.mat', nfail);
if ~isempty(L)
    R = L.R;
    Ref = table(string(R.subject), string(R.cue) + "_" + string(R.valid), string(R.band), R.dwpli, R.plv, ...
        'VariableNames', {'subject','condition','band','dwpli','plv'});
    Rh = local_roi_avg(B, 'hilbert', 'average');
    nfail = nfail + local_cmp_roi(Ref, Rh, fieldnames(cfg.battery.bands), tol_val, ...
        'T1 Hilbert/average ROI-average vs conn_results.mat');
else
    fprintf('FAIL  T1: no reference loaded\n'); nfail = nfail + 1;
end

% ---- T2 ----
[Ref, nfail] = local_ref(cfg, repo, 'imag_coherence_roi_average_gamma_high.csv', nfail);
d = local_dval(B, 'ImC', 'gamma_high', 'hilbert', 'average', conds, pairs);
d = mean(d, 2, 'omitnan'); d = d(~isnan(d));
nfail = nfail + local_cmp_dz(Ref, d, tol_dz, 'T2 Hilbert/average ImC gamma_high roi_avg d_val');

% ---- T3 ----
[Ref, nfail] = local_ref(cfg, repo, 'csd_sensitivity_summary.csv', nfail);
metrics = {'dwPLI','PLV','ImC'};
levels  = {'roi_avg','F2P4'};
for mi = 1:numel(metrics)
    dp = local_dval(B, metrics{mi}, 'gamma_high', 'hilbert', 'csd', conds, pairs);
    for li = 1:numel(levels)
        if strcmp(levels{li}, 'roi_avg'), d = mean(dp, 2, 'omitnan'); else, d = dp(:, target_pair); end
        d = d(~isnan(d));
        Rr = [];
        if ~isempty(Ref), Rr = Ref(strcmp(Ref.level, levels{li}) & strcmp(Ref.metric, metrics{mi}), :); end
        nfail = nfail + local_cmp_dz(Rr, d, tol_dz, ...
            sprintf('T3 Hilbert/CSD %s gamma_high %s d_val', metrics{mi}, levels{li}));
    end
end

% ---- T4 ----
Q = P(P.band == 'gamma_high' & P.decomposition == 'hilbert' & P.reference == 'average' & P.pair == 'F2-P4', :);
dev = abs(mod(Q.mean_angle_deg, 360) - 180);       % angular distance to 180 deg
n_within = sum(dev <= phase_tol_deg); n_cells = height(Q);
ok = n_within == phase_expect_within && n_cells == phase_expect_cells;
fprintf('\n%s  T4 F2-P4 Hilbert/average gamma_high phase difference: %d of %d cells within +-%g deg of 180 deg (expected %d of %d)\n', ...
    local_pf(ok), n_within, n_cells, phase_tol_deg, phase_expect_within, phase_expect_cells);
nfail = nfail + ~ok;

% ---- C1 ----
nfail = nfail + local_check_phase(B, P);

if nfail == 0
    fprintf('\nOVERALL: PASS\n');
else
    fprintf('\nOVERALL: FAIL (%d check(s) failed)\n', nfail);
end

% ======================= local functions =======================
function T = local_read(f, textVars, numVars)
    if ~isfile(f), error('missing %s -- run s21_battery.m first', f); end
    opts = detectImportOptions(f, 'Delimiter', ',');
    opts = setvartype(opts, textVars, 'categorical');
    opts = setvartype(opts, numVars, 'double');
    T = readtable(f, opts);
end

function R = local_roi_avg(B, dec, rf)
    % ROI average over pairs, mean(...,'omitnan') as s03
    S = B(B.decomposition == dec & B.reference == rf, :);
    [G, subject, condition, band, measure] = findgroups(S.subject, S.condition, S.band, S.measure);
    value = splitapply(@(v) mean(v, 'omitnan'), S.value, G);
    R = table(subject, condition, band, measure, value);
end

function d = local_dval(B, meas, band, dec, rf, conds, pairs)
    % [subject x pair] validity contrast, as s13/s18:
    % (avatar_valid + stick_valid)/2 - (avatar_invalid + stick_invalid)/2
    S = B(B.measure == meas & B.band == band & B.decomposition == dec & B.reference == rf, :);
    subs = unique(string(S.subject));
    [~, si] = ismember(string(S.subject), subs);
    [~, ci] = ismember(string(S.condition), conds);
    [~, pk] = ismember(string(S.pair), pairs);
    ok = ci > 0 & pk > 0;
    V = nan(numel(subs), numel(conds), numel(pairs));
    V(sub2ind(size(V), si(ok), ci(ok), pk(ok))) = S.value(ok);
    d = reshape((V(:,1,:) + V(:,3,:))/2 - (V(:,2,:) + V(:,4,:))/2, numel(subs), numel(pairs));
end

function [T, nfail] = local_ref(cfg, repo, name, nfail)
    T = [];
    f = fullfile(cfg.battery.check_dir, name);
    fprintf('\nReference: %s\n', f);
    if ~isfile(f)
        fprintf('FAIL  reference file not found\n'); nfail = nfail + 1; return;
    end
    d = dir(f);
    fprintf('  last modified: %s\n', d.date);
    alt = fullfile(repo, 'results', name);
    if isfile(alt) && ~local_same_path(f, alt)
        if local_same_content(f, alt)
            fprintf('PASS  identical to %s\n', alt);
        else
            fprintf('FAIL  contents differ from %s\n', alt); nfail = nfail + 1;
        end
    end
    [~,~,ext] = fileparts(name);
    if strcmpi(ext, '.mat')
        T = load(f);
    else
        T = readtable(f, 'TextType', 'char');
    end
end

function same = local_same_path(a, b)
    same = strcmpi(char(java.io.File(a).getCanonicalPath()), char(java.io.File(b).getCanonicalPath()));
end

function same = local_same_content(a, b)
    [~,~,ext] = fileparts(a);
    if strcmpi(ext, '.mat')
        same = isequaln(load(a), load(b));
    else
        crlf = sprintf('\r\n');
        same = strcmp(strrep(fileread(a), crlf, newline), strrep(fileread(b), crlf, newline));
    end
end

function nf = local_cmp_roi(Ref, Bat, bands, tol, label)
    % Ref: subject, condition, band, dwpli, plv (one row per subject x band x cell)
    nf = 0;
    Ref = Ref(ismember(Ref.band, bands), :);
    kr = Ref.subject + "|" + Ref.condition + "|" + Ref.band;
    meas = {'dwpli','dwPLI'; 'plv','PLV'};
    for m = 1:size(meas,1)
        Bm = Bat(Bat.measure == meas{m,2} & ismember(string(Bat.band), bands), :);
        kb = string(Bm.subject) + "|" + string(Bm.condition) + "|" + string(Bm.band);
        [tf, loc] = ismember(kr, kb);
        rv = Ref.(meas{m,1})(tf); bv = Bm.value(loc(tf));
        d = abs(rv - bv); d(isnan(rv) & isnan(bv)) = 0;
        n_missing = sum(~tf); n_extra = sum(~ismember(kb, kr));
        maxd = max([d; 0]);
        ok = n_missing == 0 && n_extra == 0 && all(d <= tol);
        fprintf('%s  %s, %s: %d rows compared, %d missing in battery, %d extra in battery, max |diff| = %.3g\n', ...
            local_pf(ok), label, meas{m,2}, sum(tf), n_missing, n_extra, maxd);
        nf = nf + ~ok;
    end
end

function nf = local_cmp_dz(Ref, d, tol, label)
    % dz = mean(d)/std(d), as s13/s18; PASS needs the same n and |diff dz| <= tol
    n = numel(d); dz = mean(d)/std(d);
    if isempty(Ref) || height(Ref) ~= 1
        fprintf('FAIL  %s: reference row not found (battery n=%d dz=%.4f)\n', label, n, dz);
        nf = 1; return;
    end
    ok = n == Ref.n && abs(dz - Ref.dz) <= tol;
    fprintf('%s  %s: battery n=%d dz=%.4f | reference n=%d dz=%.4f | |diff| = %.4f\n', ...
        local_pf(ok), label, n, dz, Ref.n, Ref.dz, abs(dz - Ref.dz));
    nf = ~ok;
end

function nf = local_check_phase(B, P)
    kB = unique(string(B.subject) + "|" + string(B.condition) + "|" + string(B.band) + "|" + ...
        string(B.decomposition) + "|" + string(B.reference) + "|" + string(B.pair));
    kP = string(P.subject) + "|" + string(P.condition) + "|" + string(P.band) + "|" + ...
        string(P.decomposition) + "|" + string(P.reference) + "|" + string(P.pair);
    n_dup = numel(kP) - numel(unique(kP));
    n_diff = numel(setxor(kB, kP));
    r = P.resultant_len;
    n_bad_r = sum(~isnan(r) & (r < 0 | r > 1 + 1e-12));
    ok = n_dup == 0 && n_diff == 0 && n_bad_r == 0;
    fprintf('\n%s  C1 battery_phasediff.csv vs battery_persubject.csv cells: %d rows, %d duplicate, %d unmatched, %d resultant_len outside [0,1]\n', ...
        local_pf(ok), numel(kP), n_dup, n_diff, n_bad_r);
    nf = ~ok;
end

function s = local_pf(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
