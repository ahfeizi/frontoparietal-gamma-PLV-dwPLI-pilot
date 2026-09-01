function Y = util_bandpass(X, srate, band, order_factor)
% Zero-phase FIR band-pass filter.
%   X    : [nchan x nsamples]
%   band : [low high] in Hz
% Returns Y same size as X. FIR (windowed-sinc) + filtfilt -> zero phase,
% the same family EEGLAB's pop_eegfiltnew uses. Requires Signal Processing Toolbox.
ord = round(order_factor * srate / band(1));
if mod(ord,2) ~= 0, ord = ord + 1; end          % even order for band-pass
bcoef = fir1(ord, band / (srate/2), 'bandpass');
Y = filtfilt(bcoef, 1, X.').';
end
