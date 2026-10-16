function [plv_pairs, dwpli_pairs, imc_pairs, mvec_pairs] = util_conn_pairs(seg, ia, ib)
% Pairwise PLV, dwPLI, imaginary coherency (ImC) and phase-difference
% resultant vector for every ia x ib channel pair of a set of analytic-signal
% segments. Used by s21_battery.m for BOTH the Hilbert and the wavelet path,
% so the estimators are identical across decompositions.
%   seg : [nchan x wlen x kept] complex analytic segments (one per trial)
%   ia  : row-channel indices into seg (frontal)
%   ib  : column-channel indices into seg (parietal)
% Pair order is p = (i-1)*numel(ib) + j, the same order as s03_connectivity.m
% and the pair labels of s13/s18.
% Returns 1 x np vectors; mvec_pairs is complex: mean(exp(1i*dphi)) pooled
% over time samples x trials.
%
% PLV / dwPLI lines are copied verbatim from s03_connectivity.m (dwPLI and PLV
% across trials per time sample, then averaged over the window).
% ImC and the phase-difference vector are copied from
% s13_imaginary_coherence.m (ImC = imag(Sxy / sqrt(Sxx*Syy)), cross- and
% auto-spectra averaged over time samples AND trials; signed).
wlen = size(seg,2); kept = size(seg,3);
np = numel(ia) * numel(ib);
dwpli_pairs = zeros(1,np); plv_pairs = zeros(1,np); imc_pairs = zeros(1,np);
mvec_pairs = complex(zeros(1,np)); p = 0;
for i = 1:numel(ia)
    zi = reshape(seg(ia(i),:,:), wlen, kept);
    for j = 1:numel(ib)
        p = p + 1;
        zj = reshape(seg(ib(j),:,:), wlen, kept);
        % ---- PLV / dwPLI (verbatim from s03_connectivity.m) ----
        X = zi .* conj(zj); im = imag(X);
        s_im = sum(im,2); s_abs = sum(abs(im),2); s_sq = sum(im.^2,2);
        dwpli_t = (s_im.^2 - s_sq) ./ (s_abs.^2 - s_sq);
        plv_t   = abs(mean(X ./ abs(X), 2));
        dwpli_pairs(p) = mean(dwpli_t, 'omitnan');
        plv_pairs(p)   = mean(plv_t, 'omitnan');
        % ---- ImC (from s13_imaginary_coherence.m) ----
        Sxy = mean(mean(zi .* conj(zj)));
        Sxx = mean(mean(abs(zi).^2));
        Syy = mean(mean(abs(zj).^2));
        coh = Sxy / sqrt(Sxx*Syy);
        imc_pairs(p) = imag(coh);
        % ---- phase-difference concentration (from s13) ----
        if nargout > 3
            dphi = angle(zi) - angle(zj);
            dphi = dphi(:);
            mvec_pairs(p) = mean(exp(1i*dphi));
        end
    end
end
end
