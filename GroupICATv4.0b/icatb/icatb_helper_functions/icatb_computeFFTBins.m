function [fftPlots, fftBins] = icatb_computeFFTBins(sesInfo, varargin)
%% Compute FFT bins
%
% Inputs:
% 1. sesInfo - sesInfo data structure
% 2. varargin - Arguments passed in pairs.
%   a. average_runs - Average runs
%   b. npointfft - N point FFT
%   c. num_bins - Number of bins
%   d. datasetno - Dataset no.
%   e. compno - Component numbers
%
% Outputs:
% 1. fftPlots - FFT of dimensions number of subjects by components by points
% 2. fftBins - FFT bins of dimensions number of subjects by components by points
%

if (ischar(sesInfo))
    load(sesInfo);
end

averageRuns = 0;
nPointFFT = 300;
numBins = 6;
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;
compno = (1:sesInfo.numComp);
icaOutputFiles = sesInfo.icaOutputFiles;

for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'average_runs'))
        averageRuns = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'npointfft'))
        nPointFFT = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'num_bins'))
        numBins = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'datasetno'))
        datasetno = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'compno'))
        compno = varargin{i + 1};
    end
end

if (~exist('datasetno', 'var'))
    if (averageRuns)
        datasetno = (1:numOfSub);
    else
        datasetno = (1:numOfSub*numOfSess);
    end
end

if (averageRuns)
    if (max(datasetno) > numOfSub)
        error('No. of subjects passed exceeds the no. of subjects in the analysis');
    end
else
    if (max(datasetno) > numOfSub*numOfSess)
        error('No. of datasets passed exceeds the no. of datasets in the analysis');
    end
end

nPoints = ceil(nPointFFT / 2);

disp('Computing fft of subject time courses ...');

subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', icaOutputFiles, 'numOfSub', numOfSub, 'numOfSess', numOfSess);

% fft and fft bins
fftPlots = zeros(length(datasetno), length(compno), nPoints);
fftBins = zeros(length(datasetno), length(compno), numBins);
datasetno = datasetno(:)';
compno = compno(:)';

tp = min(sesInfo.diffTimePoints);

countD = 0;
%% Loop over subjects
for nD = datasetno
    countD = countD + 1;
    if (averageRuns)
        tc = icatb_loadComp(sesInfo, compno, 'subjects', nD, 'sessions', (1:numOfSess), 'vars_to_load', 'tc', 'subject_ica_files', subjectICAFiles, 'average_runs', averageRuns);
    else
        subjectNumber = ceil(nD / numOfSess);
        sessNumber = mod(nD - 1, numOfSess) + 1;
        tc = icatb_loadComp(sesInfo, compno, 'subjects', subjectNumber, 'sessions', sessNumber, 'vars_to_load', 'tc', 'subject_ica_files', subjectICAFiles);
        tc = tc(1:tp, :);
    end
    
    %% Loop over components
    for nComp = 1:length(compno)
        tmp = abs(fft(icatb_normdet(tc(:, nComp), 1), nPointFFT));
        tmp = tmp(1:nPoints);
        binSize = length(tmp) / numBins;
        fftPlots(countD, nComp, :) = tmp;
        %% Loop over bins
        for nBin = 1:numBins
            % Take the average for the frequencies within each bin
            fftBins(countD, nComp, nBin) = mean(tmp((1:binSize) + binSize*(nBin - 1)));
        end
        %% End loop over bins
    end
    %% End loop over components
    
    clear tc;
    
end
%% End loop over subjects

disp('Done computing fft of subject time courses');