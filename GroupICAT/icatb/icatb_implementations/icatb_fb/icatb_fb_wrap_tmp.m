function [wt_mult, wt_mult_winidx, tmp_tc_qam] = ...
    icatb_fb_wrap_tmp(tmp_tc, winsize, wintype, modu_fs, Tr)
% Wrapper for calculate_SSB_SWPC using calculate_SSBwt interface.
%
% Inputs:
%   tmp_tc   - Time x Components
%   winsize  - Window size (odd positive integer)
%   wintype  - 'rect', 'gauss', or 'tukey'
%   modu_fs  - SSB modulation frequency (Hz)
%   Tr       - Sampling interval (seconds)
%
% Outputs:
%   wt_mult        - Windowed connectivity (Time x Pairs)
%   wt_mult_winidx - Window center indices
%   tmp_tc_qam     - SSB-modulated time courses (Time x Components)
%
% NOTE: Unlike calculate_SSBwt, wt_mult contains weighted
% Pearson correlations, not center-sample z-score products.

b_plotflag = 0;

% Call SSB sliding-window Pearson correlation
[wt_mult, wt_mult_winidx] = ...
    icatb_calculate_SSB_SWPC(tmp_tc, winsize, Tr, modu_fs, wintype, b_plotflag);

% Calculate SSB-modulated time courses for third output
if nargout > 2
    tt = (0:size(tmp_tc,1)-1)' * Tr;
    modulation_tc = exp(1i*2*pi*modu_fs*tt);
    tmp_tc_qam = real(hilbert(tmp_tc) .* modulation_tc);
end

end