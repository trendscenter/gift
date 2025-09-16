function [tValues, pValues] = icatb_compare_frequency_bins(param_file, group1, group2, compNumber)
%% Determine PSD and perform group differences
%
% Inputs:
%
% 1. param_file - ICA Parameter file
% 2. group1 - Group 1 subject numbers
% 3. group2 - Group 2 subject numbers
%
% Outputs:
%
% 1. tValues - T values of size components by (NPOINT_FFT_GROUP_COMPARISON/2).
% 2. pValues - p values of size components by (NPOINT_FFT_GROUP_COMPARISON/2).
%

%% Modality type
modalityType = icatb_get_modality;

if ~strcmpi(modalityType, 'fmri')
    error('Group comparison utility works only with the GIFT');
end

%% Load  defaults
icatb_defaults;
global PARAMETER_INFO_MAT_FILE;

% FFT defaults for group comparison
global NPOINT_FFT_GROUP_COMPARISON;
global NUM_BINS_GROUP_COMPARISON;
global DEFAULT_TR_SPECTRAL_GROUP_COMPARE;

if isempty(NPOINT_FFT_GROUP_COMPARISON)
    NPOINT_FFT_GROUP_COMPARISON = 300;
end

if isempty(NUM_BINS_GROUP_COMPARISON)
    NUM_BINS_GROUP_COMPARISON = 6;
end

if isempty(DEFAULT_TR_SPECTRAL_GROUP_COMPARE)
    DEFAULT_TR_SPECTRAL_GROUP_COMPARE = 2;
end

%% DEFAULT TR
DEFAULT_TR = DEFAULT_TR_SPECTRAL_GROUP_COMPARE;

%% N point FFT
nPointFFT = NPOINT_FFT_GROUP_COMPARISON;

%% Number of frequency bins
numBins = NUM_BINS_GROUP_COMPARISON;

%% Check bins
checkBin = mod(nPointFFT, 2*numBins);
if (checkBin ~= 0)
    error('Error:CheckBin', 'NPOINT_FFT_GROUP_COMPARISON (%d) must be exactly divisible by 2*NUM_BINS_GROUP_COMPARISON (%d).\nPlease check icatb_defaults.m file', ...
        nPointFFT, 2*numBins);
end

%% load parameters file
filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

if ~exist('param_file', 'var')
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select a valid parameter file', 'filter', filterP);
end

drawnow;

[outputDir, param_file, extn] = fileparts(param_file);

if isempty(outputDir)
    outputDir = pwd;
end

param_file = fullfile(outputDir, [param_file, extn]);

load(param_file);

if ~exist('sesInfo', 'var')
    error('Error:sesInfo', 'Selected file %s is not a valid ICA parameter file', param_file);
end

if (sesInfo.numOfSub*sesInfo.numOfSess == 1)
    error('Spectral group compare utility works with multiple datasets only');
end

sesInfo.outputDir = outputDir;

%% Get information from parameter file
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;
averageRuns = 0;

if ((numOfSub > 1) && (numOfSess > 1))
    averageRuns = icatb_questionDialog('title', 'Average across runs?', 'textbody', 'Do you want to average subject timecourses across runs?');
end

%% Select groups
group1Name = 'Group 1';
group2Name = 'Group 2';

if (averageRuns || (numOfSess == 1))
    subjectString = [repmat('Subject ', numOfSub, 1), num2str((1:numOfSub)')];
else
    if ((numOfSub > 1) && (numOfSess > 1))
        subjectString = cell(numOfSub*numOfSess, 1);
        count = 0;
        for nSub = 1:numOfSub
            for nSess = 1:numOfSess
                count = count + 1;
                subjectString{count} = ['Subject ', num2str(nSub), ' Session ', num2str(nSess)];
            end
        end
        subjectString = char(subjectString);
    else
        subjectString = [repmat('Session ', numOfSess, 1), num2str((1:numOfSess)')];
    end
end

if ~exist('group1', 'var')
    % Group Name and datasets
    [group1Name, group1] = icatb_select_groups_gui(subjectString, 'Group 1', 'select_subjects_psd');
end

if isempty(group1)
    error('Group 1 subjects are not selected');
end

if (size(subjectString, 1) < max(group1))
    error('Error:Group1', 'Maximum of group1 (%d) exceeds the maximum of subjects (%d)\n',  max(group1), size(subjectString, 1));
end

if ~exist('group2', 'var')
    % Group2 name and datasets
    [group2Name, group2] = icatb_select_groups_gui(subjectString, 'Group 2', 'select_subjects_psd');
end

if isempty(group2)
    error('Group 2 subjects are not selected');
end

if (size(subjectString, 1) < max(group2))
    error('Error:Group2', 'Maximum of group2 (%d) exceeds the maximum of subjects (%d)\n',  max(group2), size(subjectString, 1));
end

group1 = group1(:)';
group2 = group2(:)';

if (~exist('compNumber', 'var'))
    compNumber = (1:sesInfo.numComp);
end

compNumber = compNumber(:)';

drawnow;

%% Print some information
disp(sprintf('N point FFT used for group comparison is %d', nPointFFT));
disp(sprintf('Number of frequency bins used is %d', numBins));
disp(sprintf('Selected subjects for group 1 are %s ', num2str(group1(:)')));
disp(sprintf('Selected subjects for group 2 are %s ', num2str(group2(:)')));

%% End for selecting groups

%nPointFFT = timePoints(1);
nPoints = ceil(nPointFFT / 2);

fprintf('\n');

%% Compute FFT bins
[fftPlots, fftBins] = icatb_computeFFTBins(sesInfo, 'average_runs', averageRuns, 'npointfft', nPointFFT, 'num_bins', numBins, 'datasetno', [group1, group2], 'compno', compNumber);

%% Compute mean power of each group
mean_power_group1 = squeeze(mean(fftPlots(1:length(group1), :, :)));
mean_power_group2 = squeeze(mean(fftPlots(length(group1) + 1:length(group1) + length(group2), :, :)));

%clear fftPlots;

fprintf('\n');

%% Create two matrices for between-group comparisons.

disp(['Comparing frequency bins of groups (', group1Name, ' & ', group2Name, ') ...']);

group1FFTBins = fftBins(1:length(group1), :, :);
group2FFTBins = fftBins(length(group1) + 1:length(group1) + length(group2), :, :);

binSize = nPoints/numBins;

%% Total number of subjects
numSampleA = length(group1);
numSampleB =  length(group2);
df = (numSampleA + numSampleB - 2);

%% Initialise tValues and pValues
tValues = zeros(size(group1FFTBins, 2), size(group1FFTBins, 3));
pValues = zeros(size(group1FFTBins, 2), size(group1FFTBins, 3));

%% Loop over components
for nComp = 1:size(group1FFTBins, 2)
    
    %% Sample A and B
    sampleA = squeeze(group1FFTBins(:, nComp,:));
    sampleB = squeeze(group2FFTBins(:, nComp,:));
    
    a = (numSampleA - 1).*(std(sampleA).^2);
    b = (numSampleB - 1).*(std(sampleB).^ 2);
    
    %% Numerator
    num = (sqrt(numSampleA*numSampleB)).*(squeeze(mean(sampleA)) - squeeze(mean(sampleB)));
    
    %% Denominator
    denom = sqrt((a + b).*(numSampleA + numSampleB) ./ (numSampleA + numSampleB - 2));
    
    clear a b;
    
    %% T and p values (Two tailed t-test)
    tValues(nComp, :) = num./denom;
    
    pValues(nComp, :) = icatb_get_pvalue(tValues(nComp, :), df, 0);
    
    clear num denom;
    
end
%% End loop over components

disp('Done comparing frequency bins');

fprintf('\n');

%% XTICK Label
TR = DEFAULT_TR;
freq = (1/TR/2);
xtickVal = linspace(0, freq, nPoints);
xtickVal = round(xtickVal*100)/100;
inds = ceil(linspace(binSize, nPoints, numBins));
if (inds(end) > length(xtickVal))
    inds(end) = length(xtickVal);
end

freqValues = xtickVal(inds);
freqValues = freqValues(:);

xtickLabel = cellstr(num2str(freqValues));
%xtickLabel = repmat(xtickLabel, size(tValues, 2), 1);

%% Plot Group 1 - Group 2
%Multi bar plot
plotTitle = {[group1Name, ' - ', group2Name]};

if (size(group1FFTBins, 2) > 4)
    elemPerAxes = 4;
else
    elemPerAxes = size(group1FFTBins, 2);
end

legendStr = cellstr([repmat('Comp ', size(group1FFTBins, 2), 1), num2str(compNumber')]);

figTitle = 'Power Spectrum Analysis of ICA timecourses';
icatb_multi_bar_plot(tValues', 'XTickLabel', xtickLabel, 'xlabel', 'Frequency (Hz)', 'ylabel', 'T-Values', 'title', figTitle, ...
    'axesTitle', plotTitle, 'elemPerAxes', elemPerAxes, 'legend', legendStr);

%% Save results
outFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_comparison_frequency_bins.mat']);
icatb_save(outFile, 'tValues', 'pValues', 'group1Name', 'group2Name', 'mean_power_group1', 'mean_power_group2', 'fftPlots', 'fftBins');
disp(['Results are saved in file ', outFile]);

fprintf('\n');
