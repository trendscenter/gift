function total_percent_variance  = icatb_percent_variance(param_file)
%% Use regression to get explained percent variance in the data by the model
% where components are treated as model.
%
% Inputs:
% 1. param_file - Parameter file.
%
% Outputs:
% 1. percent_variance - Percent variance

global PARAMETER_INFO_MAT_FILE;

%% Load parameter and validate parameter file
filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];

if ~exist('param_file', 'var')
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);
end

if isempty(param_file)
    error('Error:ParameterFile', 'Parameter file is not selected\n');
end

drawnow;

load(param_file);

outputDir = fileparts(param_file);

if isempty(outputDir)
    outputDir = pwd;
end

cd(outputDir);

if ~exist('sesInfo', 'var')
    error('Error:ParameterFile', 'Selected file %s is not a valid parameter file\n', param_file);
end

if (~sesInfo.isInitialized)
    error('Please run the analysis in order to compute percent variance');
end

[modalityType, dataTitle, compSetFields] = icatb_get_modality;

%% Get required vars from parameter file
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;

%% Initialise variables
count = 0;
ss_residual = 0;
ss_total= 0;

fprintf('\n');
disp('Computing percent variance explained by the components in the data ...');
fprintf('\n');
disp('This part will take time depending on the number of subjects and components. Please wait ...');
fprintf('\n');

%% Loop over subjects
for nSub = 1:numOfSub
    %% Loop over sessions
    for nSess = 1:numOfSess
        count = count + 1;

        % Load back-reconstruction MAT file
        subFile = [sesInfo.back_reconstruction_mat_file, num2str(count), '.mat'];
        subFile = fullfile(outputDir, subFile);
        if (exist(subFile, 'file'))
            load(subFile);
        else
            sesInfo.conserve_disk_space = 1;
            compSet = getBackReconSet(sesInfo, count);
        end

        % Mixing matrix
        tc = getfield(compSet(1), compSetFields{2});
        clear compSet

        %% Load data and apply mask
        files = char(sesInfo.inputFiles(count).name);

        data = icatb_read_data(files, [], sesInfo.mask_ind);

        % Compute sum of squares of residual and total sum of squares
        [temp_ss_residual, temp_ss_total]  = icatb_compute_percent_variance(data, tc);

        clear tc data;

        % Add the current sum of squares residual to total
        ss_residual = ss_residual + temp_ss_residual;
        ss_total = ss_total + temp_ss_total;

    end
    %% End loop over sessions
end
%% End loop over subjects

%% Percent variance in the data
total_percent_variance = 100*(1 - (ss_residual/ss_total));

disp(['The total percent variance explained by components in the data is ', num2str(total_percent_variance)]);

fprintf('\n');

disp('Done');
fprintf('\n');

function compSet = getBackReconSet(sesInfo, nDataSet)
%% Get back-reconstructed component of the subject
%

sesInfo.dataSetNo = nDataSet;
[sesInfo, compSet] = icatb_backReconstruct(sesInfo);