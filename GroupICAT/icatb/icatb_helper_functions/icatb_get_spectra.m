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

icatb_defaults;
global TIMECOURSE_POSTPROCESS;

if (length(data) == numel(data))
    data = data(:)';
end

if (~exist('TR', 'var'))
    TR = 1;
end

NPointFFT = 300;
NBins = 6;
option = 1;

try
    option = TIMECOURSE_POSTPROCESS.spectra.option;
catch
end

try
    NPointFFT = TIMECOURSE_POSTPROCESS.spectra.NPointFFT;
    NBins = TIMECOURSE_POSTPROCESS.spectra.NBins;
catch
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
    if (option == 1)
        [tmp, f] = icatb_mtspectrumc(data(nSub, :), params);
    else
        %[tmp, f] = computeFFT(data(nSub, :), NPointFFT, NBins);
        [tmp, f] = computeFFT(data(nSub, :), TR, NPointFFT, NBins);
    end
    if (nSub == 1)
        spectra_tc = zeros(size(data, 1), length(tmp));
    end
    spectra_tc(nSub, :) = tmp(:)';
end


function [fftBins, freq] = computeFFT(tc, TR, NPointFFT, NBins)
% Compute spectra over bins
%

nPoints = ceil(NPointFFT / 2);

tmp = abs(fft(icatb_normdet(tc, 1), NPointFFT));
tmp = tmp(1:nPoints);
binSize = length(tmp) / NBins;

fftBins = zeros(1, NBins);
%% Loop over bins
for nBin = 1:NBins
    % Take the average for the frequencies within each bin
    fftBins(nBin) = mean(tmp((1:binSize) + binSize*(nBin - 1)));
end

freq = linspace(0, (1/TR/2), NBins);