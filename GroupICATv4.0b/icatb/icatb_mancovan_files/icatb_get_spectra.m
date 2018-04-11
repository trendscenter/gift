function [spectra_tc, f] = icatb_get_spectra(data, TR, params)
%% Get multi-taper spectral estimation
%
% Inputs:
% 1. data - 2D array of dimensions subjects by timepoints
% 2. TR - TR of the experiment in seconds
%
% Outputs:
% 1. spectra_tc - Spectra
% 2. f - Freq in Hz
%

if (length(data) == numel(data))
    data = data(:)';
end

if (~exist('TR', 'var'))
    TR = 1;
end

defaultParams = struct('tapers', [3, 5], 'Fs', (1/TR), 'fpass', [0, 1/(2*TR)]);

if (~exist('params', 'var'))
    params = defaultParams;
end

if (~isfield(params, 'tapers'))
    params.tapers = defaultParams.tapers;
end

if (~isfield(params, 'Fs'))
    params.Fs = defaultParams.Fs;
end

if (~isfield(params, 'fpass'))
    params.fpass = defaultParams.fpass;
end

for nSub = 1:size(data, 1)
    [tmp, f] = icatb_mtspectrumc(data(nSub, :), params);
    if (nSub == 1)
        spectra_tc = zeros(size(data, 1), length(tmp));
    end
    spectra_tc(nSub, :) = tmp(:)';
end