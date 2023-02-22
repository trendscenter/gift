function [dynamicrange, fALFF] = icatb_get_spec_stats(s, f, freq_limits)
% Get spectral stats
%

if (~exist('freq_limits', 'var'))
    freq_limits = [0.1, 0.15];
end

lower_limit = min(freq_limits);
upper_limit = max(freq_limits);


[mv, ind]=max(s);
dynamicrange = max(s)-min(s(ind:end));
LF = find(f < lower_limit);
HF = find(f > upper_limit);
fALFF = trapz(f(LF),s(LF))/trapz(f(HF),s(HF));