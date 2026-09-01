% s04_stats.m  --  Group-level within-subject contrasts, Pilot 2.
% Unlike Pilot 1's single two-level factor (load), this dataset has a 2x2
% within-subject design: cue (avatar/stick) x validity (valid/invalid).
% For each band and each metric (dwpli, plv, frontal_power, parietal_power),
% tests THREE contrasts via non-parametric sign-flip permutation:
%   cue_main   : mean(avatar) - mean(stick)               [averaged over validity]
%   valid_main : mean(valid)  - mean(invalid)              [averaged over cue]
%   interaction: (avatar_valid - avatar_invalid) - (stick_valid - stick_invalid)
% Reports mean per condition/level, paired diff, Cohen's dz, t, and
% permutation p for each contrast.

clear; clc;
cfg = config();
load(fullfile(cfg.deriv_dir,'conn_results.mat'), 'R');
rng(1);

bands   = fieldnames(cfg.bands);
metrics = {'dwpli','plv','frontal_power','parietal_power'};
subs    = unique(R.subject);
results = struct();

% helper: pull the 4-cell value for one subject/band/metric/cue/valid, or NaN
getval = @(subj,bn,mt,cueLab,valLab) local_getval(R,subj,bn,mt,cueLab,valLab);

fprintf('\n band       metric          contrast     n     A        B      dz      t       p\n');
fprintf('----------------------------------------------------------------------------------\n');

for b = 1:numel(bands)
  bn = bands{b};
  for m = 1:numel(metrics)
    mt = metrics{m};
    use_log = cfg.log_power && ismember(mt, {'frontal_power','parietal_power'});

    % pull all 4 cells for every subject
    av = nan(numel(subs),1); ai = nan(numel(subs),1);
    sv = nan(numel(subs),1); si = nan(numel(subs),1);
    for i = 1:numel(subs)
        av(i) = getval(subs{i}, bn, mt, 'avatar', 'valid');
        ai(i) = getval(subs{i}, bn, mt, 'avatar', 'invalid');
        sv(i) = getval(subs{i}, bn, mt, 'stick',  'valid');
        si(i) = getval(subs{i}, bn, mt, 'stick',  'invalid');
    end
    if use_log
        av = log(av); ai = log(ai); sv = log(sv); si = log(si);
    end

    contrasts = struct();
    contrasts.cue_main    = struct('A', (av+ai)/2, 'B', (sv+si)/2, 'labA','avatar','labB','stick');
    contrasts.valid_main  = struct('A', (av+sv)/2, 'B', (ai+si)/2, 'labA','valid','labB','invalid');
    contrasts.interaction = struct('A', av-ai, 'B', sv-si, 'labA','(avat: v-i)','labB','(stick: v-i)');

    cnames = fieldnames(contrasts);
    for c = 1:numel(cnames)
        cn = cnames{c};
        Avals = contrasts.(cn).A; Bvals = contrasts.(cn).B;
        ok = ~isnan(Avals) & ~isnan(Bvals);
        d = Avals(ok) - Bvals(ok); n = sum(ok);
        if n < 3
            fprintf('%6s %-14s %-12s  n=%d too few\n', bn, mt, cn, n);
            continue;
        end

        tobs = mean(d) / (std(d)/sqrt(n));
        np = cfg.stat.n_perm; tnull = zeros(np,1);
        for pp = 1:np
            dp = d .* sign(rand(n,1)-0.5);
            tnull(pp) = mean(dp) / (std(dp)/sqrt(n));
        end
        switch cfg.stat.tail
          case 0,  pval = mean(abs(tnull) >= abs(tobs));
          case 1,  pval = mean(tnull >= tobs);
          case -1, pval = mean(tnull <= tobs);
        end

        results.(bn).(mt).(cn) = struct('n',n, 'meanA',mean(Avals(ok)), 'meanB',mean(Bvals(ok)), ...
            'mean_diff',mean(d), 'cohen_dz',mean(d)/std(d), 't',tobs, 'p',pval, ...
            'labA',contrasts.(cn).labA, 'labB',contrasts.(cn).labB);

        fprintf('%6s %-14s %-12s %3d  %7.4f  %7.4f  %5.2f  %6.2f  %.4f\n', ...
            bn, mt, cn, n, mean(Avals(ok)), mean(Bvals(ok)), mean(d)/std(d), tobs, pval);
    end
  end
end

save(fullfile(cfg.deriv_dir,'stats_results.mat'), 'results');
fprintf('\nSaved stats_results.mat\n');

% ======================= local functions =======================
function v = local_getval(R, subj, bn, mt, cueLab, valLab)
    r = R(strcmp(R.subject,subj) & strcmp(R.band,bn) & ...
          strcmp(R.cue,cueLab) & strcmp(R.valid,valLab), :);
    if isempty(r), v = NaN; else, v = r.(mt)(1); end
end
