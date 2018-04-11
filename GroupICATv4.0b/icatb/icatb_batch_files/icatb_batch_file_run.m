function icatb_batch_file_run(inputFiles)
%% Batch file for running group ICA
%
% Inputs:
% inputFiles - location of the inputFiles (fullfile path incase file is not
% on MATLAB search path)
%
%

if (~exist('inputFiles', 'var'))
    inputFiles = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'multiple', 'filter', '*.m', 'title', 'Select input batch file/files ...');
    drawnow;
end

if (isempty(inputFiles))
    error('Input M file is not selected for batch analysis');
end

%% Turn off finite warning
warning('off', 'MATLAB:FINITE:obsoleteFunction');

%% Turn off nifti class warning when loading SPM mat files
warning('off', 'MATLAB:unknownElementsNowStruc');

warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');

warning('off', 'MATLAB:pfileOlderThanMfile');

%% Delete previous figures
%icatb_delete_gui({'groupica', 'eegift', 'gift'});

%% Check version and run this to display pushbuttons with right background
% on Matlab version 14 and later
try
    % add this statement to fix the button color
    feature('JavaFigures', 0);
catch
end

%% Open Matlab server
icatb_openMatlabServer;

inputFiles = cellstr(inputFiles);

inputFiles = formFullPaths(inputFiles);

for nFile = 1:length(inputFiles)
    
    %% Set up ICA
    param_file = icatb_setup_analysis(inputFiles{nFile});
    
    load(param_file);
    
    %% Run Analysis (All steps)
    icatb_runAnalysis(sesInfo, 1);
    
    clear sesInfo;
    
    display_results = 0;
    try
        display_results = sesInfo.display_results;
    catch
    end
    
    if ((display_results ~= 0) && ~strcmpi(sesInfo.modality, 'eeg'))
        results.formatName = display_results;
        icatb_report_generator(param_file, results);
    end
    
end

function inputFiles = formFullPaths(inputFiles)
%% Form full paths

oldDir = pwd;

for nFile = 1:length(inputFiles)
    
    cF = inputFiles{nFile};
    [p, fN, extn] = fileparts(deblank(cF));
    
    if (isempty(p))
        p = fileparts(which(cF));
        if (isempty(p))
            error('Error:InputFile', 'File %s doesn''t exist\n', cF);
        end
    end
    
    cd(p);
    
    inputFiles{nFile} = fullfile(pwd, [fN, extn]);
    
    cd(oldDir);
    
end
