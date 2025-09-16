function msk_inds = icatb_window_mask(Ashift2, stepNo, window_alpha, wsize)
%% Window mask
%
% Remove values outside the window
%

icatb_defaults;
global DFNC_DEFAULTS;

nT = length(Ashift2);

mask_windows = 1;
try
    mask_windows = DFNC_DEFAULTS.mask_windows;
catch
end

zero_val = 1e-4;

if (mask_windows)
    tmp_mask = zeros(length(Ashift2),1);
    minTp = max([1, stepNo - round((window_alpha/2)*wsize)]);
    maxTp = min([nT, stepNo + round((window_alpha/2)*wsize)]);
    tmp_mask(minTp:maxTp)=1;
    msk_inds = find(tmp_mask==1);
    msk_inds(find(Ashift2(msk_inds) <= zero_val))=[];
else
    msk_inds = (1:nT);
end