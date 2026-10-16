% s24_imc_diagnostic_check.m -- Regression check of s23_imc_diagnostic.m.
% imc_joint and dwpli_trial in imc_diag_persubject.csv must equal the ImC
% and dwPLI rows of battery_persubject.csv (band gamma_high) for every
% subject x condition x decomposition x reference x pair, to 1e-6.
% PASS needs every diagnostic row to have a battery row and vice versa,
% and max |diff| <= 1e-6 (NaN in both counts as equal).
% Both files are read from cfg.battery.check_dir (default: <repo>/results/;
% set it to cfg.deriv_dir to check a fresh s23 run).

clear; clc;
repo = fileparts(fileparts(mfilename('fullpath')));   % repository root (this script is in battery/)
addpath(repo, fullfile(repo, 'battery'));             % config.m, util_bandpass.m, util_chanidx.m; util_conn_pairs.m
cfg = config();
tol_val = 1e-6;
nfail = 0;

fprintf('Reading outputs from %s ...\n', cfg.battery.check_dir);
D = local_read(fullfile(cfg.battery.check_dir,'imc_diag_persubject.csv'), ...
    {'subject','condition','decomposition','reference','pair','measure'}, {'value'});
B = local_read(fullfile(cfg.battery.check_dir,'battery_persubject.csv'), ...
    {'subject','condition','band','measure','decomposition','reference','pair'}, {'value'});
B = B(B.band == 'gamma_high', :);

pairs = {'imc_joint','ImC'; 'dwpli_trial','dwPLI'};
for k = 1:size(pairs,1)
    Dk = D(D.measure == pairs{k,1}, :);
    Bk = B(B.measure == pairs{k,2}, :);
    kd = string(Dk.subject) + "|" + string(Dk.condition) + "|" + string(Dk.decomposition) + "|" + ...
        string(Dk.reference) + "|" + string(Dk.pair);
    kb = string(Bk.subject) + "|" + string(Bk.condition) + "|" + string(Bk.decomposition) + "|" + ...
        string(Bk.reference) + "|" + string(Bk.pair);
    [tf, loc] = ismember(kd, kb);
    dv = Dk.value(tf); bv = Bk.value(loc(tf));
    d = abs(dv - bv); d(isnan(dv) & isnan(bv)) = 0;
    n_missing = sum(~tf); n_extra = sum(~ismember(kb, kd));
    n_dup = numel(kd) - numel(unique(kd));
    maxd = max([d; 0]);
    ok = n_missing == 0 && n_extra == 0 && n_dup == 0 && all(d <= tol_val);
    fprintf('%s  %s vs battery %s (gamma_high): %d rows compared, %d not in battery, %d battery rows not in diagnostic, %d duplicate, max |diff| = %.3g\n', ...
        local_pf(ok), pairs{k,1}, pairs{k,2}, sum(tf), n_missing, n_extra, n_dup, maxd);
    nfail = nfail + ~ok;
end

if nfail == 0
    fprintf('\nOVERALL: PASS\n');
else
    fprintf('\nOVERALL: FAIL (%d check(s) failed)\n', nfail);
end

% ======================= local functions =======================
function T = local_read(f, textVars, numVars)
    if ~isfile(f), error('missing %s', f); end
    opts = detectImportOptions(f, 'Delimiter', ',');
    opts = setvartype(opts, textVars, 'categorical');
    opts = setvartype(opts, numVars, 'double');
    T = readtable(f, opts);
end

function s = local_pf(ok)
    if ok, s = 'PASS'; else, s = 'FAIL'; end
end
