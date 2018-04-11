function sesInfo = icatb_runAnalysis(sesInfo, groupICAStep)
%% Run analysis consists of the following options:
% 1. All*** - Run all the steps involved in group ICA
% 2. Initialize Parameters - Parameters are intialized after the setup
% ICA analysis
% 3. Group Data Reduction - fMRI data is reduced using the data reduction
% step
% 4. Calculate ICA - Involves calculating ICA - Several Algorithms are
% available
% 5. Back Reconstruct - Performs back reconstruction step
% 6. Calibrate Components - Calibrates the arbitrary parameters to
% percent signal changes
% 7. Group Stats - Perform stats for a group of subjects. For one subject
% one session stats are not performed.
% 8. Resume - Resumes the interrupted analysis.
%

% Inputs:
% 1. sesInfo - structure containing all parameters for analysis
% 2. groupICAStep - Step no.
%   1 - All***
%   2 - Initialize parameters
%   3 - Data Reduction
%   4 - Calculate ICA
%   5 - Back Reconstruct
%   6 - Calibrate Components
%   7 - Group Stats
%   8 - Resume


try
    % defaults
    icatb_defaults;
    
    %Screen Color Defaults
    global BG_COLOR;
    global FONT_COLOR;
    global AXES_COLOR;
    global PARAMETER_INFO_MAT_FILE;
    global ZIP_IMAGE_FILES;
    global OPEN_DISPLAY_GUI;
    global SPM_STATS_WRITE_TAL;
    global SPM_STATS_AVG_RUNS;
    
    % Return modality type
    modalityType = icatb_get_modality;
    
    open_display_gui = OPEN_DISPLAY_GUI;
    
    % load parameters file
    filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
    
    analysisType = 'batch';
    
    % load valid parameter file
    if ~exist('sesInfo', 'var')
        % show directions about run ica
        [P] = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);
        if isempty(P)
            error('Parameter file is not selected for analysis');
        end
        [pathstr, file] = fileparts(P);
        outputDir = pathstr; % output directory
        cd(pathstr);
        load(P);
        if ~exist('sesInfo', 'var')
            error('Not a valid parameter file');
        else
            disp('Parameters file succesfully loaded');
        end
        %% output directory
        sesInfo.outputDir = outputDir; % set the output directory
    else
        outputDir = sesInfo.userInput.pwd;
    end
    
    drawnow;
    
    if isfield(sesInfo.userInput, 'modality')
        if ~strcmpi(modalityType, sesInfo.userInput.modality)
            if strcmpi(sesInfo.userInput.modality, 'fmri')
                error('Use GIFT toolbox to run the analysis on MRI data.');
            elseif strcmpi(sesInfo.userInput.modality, 'fmri')
                error('Use EEGIFT toolbox to run the analysis on EEG data.');
            else
                error('Use SBM toolbox to run the analysis on sMRI data.');
            end
        end
    else
        sesInfo.userInput.modality = modalityType;
    end
    
    spm_stats_write_tal = 0;
    if (~isempty(SPM_STATS_WRITE_TAL))
        spm_stats_write_tal = SPM_STATS_WRITE_TAL;
    end
    
    spm_stats_avg_runs = 0;
    if (~isempty(SPM_STATS_AVG_RUNS))
        spm_stats_avg_runs = SPM_STATS_AVG_RUNS;
    end
    
    doSPMStats = 0;
    
    if (strcmpi(modalityType, 'fmri'))
        if (sesInfo.userInput.numOfSub > 1)
            doSPMStats = 1;
        elseif ((sesInfo.userInput.numOfSub == 1) && (sesInfo.userInput.numOfSess > 1))
            if (~spm_stats_avg_runs)
                doSPMStats = 1;
            end
        end
    end
    
    doSPMStats = doSPMStats && spm_stats_write_tal;
    
    if (doSPMStats)
        spmPath = which('spm.m');
        
        if isempty(spmPath);
            error('SPM does not exist on MATLAB path. Set SPM_STATS_WRITE_TAL to 0 if you don''t want to do SPM Stats');
        end
        
        verNum = str2num(strrep(lower(spm('ver')), 'spm', ''));
        
        if (verNum < 5)
            error('SPM stats utility works with SPM5 and higher');
        end
    end
    
    sesInfo.userInput.pwd = outputDir;
    sesInfo.outputDir = outputDir;
    
    %% Group PCA settings
    perfOptions = icatb_get_analysis_settings;
    perfType = 'user specified settings';
    if (isfield(sesInfo.userInput, 'perfType'))
        perfType = sesInfo.userInput.perfType;
    end
    
    if (isnumeric(perfType))
        perfType = perfOptions{perfType};
    end
    
    perfType = lower(perfType);
    
    ica_types = cellstr(icatb_icaAlgorithm);
    if (~ischar(sesInfo.userInput.algorithm))
        algorithmName = ica_types{sesInfo.userInput.algorithm};
    else
        algorithmName = sesInfo.userInput.algorithm;
    end
    
    if (strcmpi(algorithmName, 'moo-icar'))
        algorithmName = 'gig-ica';
    end
    
    analysisStr = {'All***', 'Initialize Parameters', 'Group Data Reduction', 'Calculate ICA/IVA', 'Back Reconstruct', 'Calibrate Components', 'Group Stats'};
    
    allSteps = {'all', 'parameter_initialization', 'group_pca', 'calculate_ica', 'back_reconstruct', 'scale_components', 'group_stats', 'resume'};
    
    % choose which steps to preform for group ica
    if ~exist('groupICAStep', 'var')
        
        disp('Opening run analysis GUI. Please wait ...');
        
        run_analysis_settings = icatb_run_analysis_settings(sesInfo);
        
        perfType = run_analysis_settings.perfType;
        stepsToRun = run_analysis_settings.stepsToRun;
        analysisType = 'GUI';
        
    else
        
        stepsToRun = groupICAStep;
        
    end
    
    if (~isnumeric(stepsToRun))
        stepsToRun = lower(cellstr(stepsToRun));
        [dd, stepsToRun] = intersect(allSteps, stepsToRun);
    end
    
    stepsToRun = sort(unique(stepsToRun));
    
    if any(stepsToRun == 1)
        stepsToRun = (2:7);
        userInput = sesInfo.userInput;
        outputDir = sesInfo.outputDir;
        if (isfield(sesInfo, 'zipContents'))
            zipContents = sesInfo.zipContents;
        end
        sesInfo = [];
        sesInfo.userInput = userInput;
        sesInfo.outputDir = outputDir;
        if (exist('zipContents', 'var'))
            sesInfo.zipContents = zipContents;
        end
        clear userInput;
    elseif any(stepsToRun == 8)
        [resume_info, sesInfo] = icatb_get_resume_info(sesInfo);
        if (isempty(resume_info))
            return;
        else
            stepsToRun = resume_info.stepsToRun;
            reductionStepNo = resume_info.groupNo;
            dataSetNo = resume_info.dataSetNo;
            val = icatb_questionDialog('title', 'Resume Analysis', 'textbody', sprintf(['Toolbox detected part/parts of the analysis to be run. Do you want to run the following step/steps?\n', repmat('\n  %s', 1, length(stepsToRun))], analysisStr{stepsToRun}));
            if (~val)
                return;
            end
            
            %% Save parameter file
            [pp, fileName] = fileparts(sesInfo.userInput.param_file);
            icatb_save(fullfile(outputDir, [fileName, '.mat']), 'sesInfo');
            clear fileName;
            
        end
        clear resume_info;
    end
    
    
    stepsToRun = stepsToRun(:)';
    
    % check the user input
    if (any(stepsToRun == 2))
        if isfield(sesInfo, 'reduction')
            sesInfo = rmfield(sesInfo, 'reduction');
        end
    end
    
    conserve_disk_space = 0;
    if (isfield(sesInfo.userInput, 'conserve_disk_space'))
        conserve_disk_space = sesInfo.userInput.conserve_disk_space;
    end
    
    %% Open parallel mode
    parallel_info.mode = 'serial';
    parallel_info.num_workers = 4;
    
    try
        parallel_info = sesInfo.userInput.parallel_info;
    catch
    end
    
    parallelMode = parallel_info.mode;
    num_workers = parallel_info.num_workers;
    
    useTemporalICA = 0;
    if (strcmpi(modalityType, 'fmri'))
        try
            useTemporalICA = strcmpi(sesInfo.userInput.group_ica_type, 'temporal');
        catch
        end
    end
    
    if (conserve_disk_space == 1)
        stepsToRun(stepsToRun == 5) = [];
        stepsToRun(stepsToRun == 7) = [];
    end
    
    if (useTemporalICA)
        stepsToRun(stepsToRun == 5) = [];
    end
    
    
    if (strcmpi(algorithmName, 'gig-ica') || strcmpi(algorithmName, 'constrained ica (spatial)'))
        % No data reduction
        stepsToRun(stepsToRun == 3) = [];
        % No back-reconstruction
        stepsToRun(stepsToRun == 5) = [];
    end
    
    
    if (isfield(sesInfo.userInput, 'modality'))
        sesInfo.modality = sesInfo.userInput.modality;
    else
        sesInfo.modality = modalityType;
    end
    
    % performing batch analysis
    output_LogFile = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, '_results.log']);
    
    % Print output to a file
    diary(output_LogFile);
    
    toolboxNames = ver;
    parallelCluster= ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));
    
    sesInfo.parallel_info = parallel_info;
    
    if (strcmpi(parallelMode, 'parallel') || strcmpi(analysisType, 'batch'))
        statusHandle = [];
        disp('Starting Analysis ');
        fprintf('\n');
    else
        % GUI
        
        waitbarTag = [sesInfo.userInput.prefix, 'waitbar'];
        waitbarCheckH = findobj('tag', waitbarTag);
        try
            delete(waitbarCheckH);
        catch
        end
        
        % Include a status bar that shows what percentage of analysis is
        % completed
        titleFig = 'System busy. Please wait...'; perCompleted = 0;
        statusHandle = waitbar(perCompleted, titleFig, 'name', [num2str(perCompleted*100), '% analysis done'], ...
            'DefaultTextColor', FONT_COLOR, 'DefaultAxesColor', AXES_COLOR, 'color', BG_COLOR, 'tag', waitbarTag);
        
        appDataName = 'gica_waitbar_app_data';
        if (isappdata(statusHandle, appDataName))
            rmappdata(statusHandle, appDataName);
        end
        
        % Number of calls per function where most of the time is spent in
        % analysis
        numberOfCalls_function = 3;
        
        unitPerCompleted = 1/(length(stepsToRun)*numberOfCalls_function);
        
        
        statusData.unitPerCompleted = unitPerCompleted;
        statusData.perCompleted = 0;
        
        setappdata(statusHandle, appDataName, statusData);
        
    end
    
    
    sesInfo.userInput.perfType = perfType;
    
    % Use tic and toc instead of cputime
    tic;
    
    if (strcmpi(parallelMode, 'parallel') && parallelCluster)
        if (~isempty(which('parpool')))
            try
                parpool(num_workers);
            catch
            end
        else
            try
                matlabpool('open', num_workers);
            catch
            end
        end
    end
    
    subjectFile = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix, 'Subject.mat']);
    
    if (~exist(subjectFile, 'file'))
        files = sesInfo.userInput.files;
        numOfSub = sesInfo.userInput.numOfSub;
        numOfSess = sesInfo.userInput.numOfSess;
        SPMFiles = sesInfo.userInput.designMatrix;
        icatb_save(subjectFile, 'files', 'numOfSub', 'numOfSess', 'SPMFiles', 'modalityType');
        clear files SPMFiles numOfSub numOfSess;
    end
    
    
    
    countStep = 0;
    for groupICAStep = stepsToRun
        
        countStep = countStep + 1;
        
        % parameter initialization
        if(groupICAStep == 1 || groupICAStep == 2)
            sesInfo = icatb_parameterInitialization(sesInfo, statusHandle);
        end
        
        if ((countStep == 1) && exist('dataSetNo', 'var'))
            sesInfo.dataSetNo = dataSetNo;
        end
        
        if (~strcmpi(modalityType, 'eeg'))
            if (isfield(sesInfo.HInfo.V(1), 'private') && ~isa(sesInfo.HInfo.V(1).private, 'icatb_nifti'))
                %if (~isa(sesInfo.HInfo.V(1).private, 'icatb_nifti'))
                [dd, sesInfo.HInfo] = icatb_returnHInfo(sesInfo.HInfo.V(1).fname);
                sesInfo.userInput.HInfo = sesInfo.HInfo;
            end
        end
        
        % Data reduction
        if(groupICAStep == 1 || groupICAStep == 3)
            
            if (exist('reductionStepNo', 'var'))
                reductionStepsToRun = (1:sesInfo.numReductionSteps);
                reductionStepsToRun(reductionStepsToRun < reductionStepNo) = [];
                sesInfo.reductionStepsToRun = reductionStepsToRun;
            end
            
            if ((strcmpi(perfType, 'maximize performance')) || (strcmpi(perfType, 'less memory usage')))
                
                % Get performance settings
                [max_mem, pcaType, pcaOpts] = icatb_get_analysis_settings(sesInfo, perfType);
                
                gpca_opts = cell(1, sesInfo.numReductionSteps);
                if (strcmpi(algorithmName, 'iva-gl') || strcmpi(algorithmName, 'iva-l') || (sesInfo.numReductionSteps > 1))
                    % use double precision with svd to handle
                    % ill-conditioned datasets
                    tmp_pca_opts = pcaOpts;
                    try
                        tmp_pca_opts.precision = 'double';
                        tmp_pca_opts.solver = 'all';
                    catch
                    end
                    gpca_opts{1} = struct('pcaType', 'svd', 'pca_opts', icatb_pca_options('svd', tmp_pca_opts, 'off'));
                else
                    gpca_opts{1} = struct('pcaType', pcaType, 'pca_opts', icatb_pca_options(pcaType, pcaOpts, 'off'));
                end
                
                if (sesInfo.numReductionSteps > 1)
                    gpca_opts(2:sesInfo.numReductionSteps) = repmat({struct('pcaType', pcaType, 'pca_opts', pcaOpts)}, 1, sesInfo.numReductionSteps - 1);
                end
                
                sesInfo.pca_opts = gpca_opts;
                
            end
            
            sesInfo = icatb_dataReduction(sesInfo, statusHandle);
            
        end
        
        % Calculate ICA
        if(groupICAStep == 1 || groupICAStep == 4)
            sesInfo = icatb_calculateICA(sesInfo, statusHandle);
        end
        
        % Back Reconstruction
        if(groupICAStep == 1 || groupICAStep == 5)
            sesInfo = icatb_backReconstruct(sesInfo, statusHandle);
        end
        
        % Calibrate Components
        if(groupICAStep == 1 || groupICAStep == 6)
            sesInfo = icatb_calibrateComponents(sesInfo, statusHandle);
        end
        
        % Group Stats
        if(groupICAStep == 1 || groupICAStep == 7)
            sesInfo = icatb_groupStats(sesInfo, statusHandle);
        end
        
    end
    
    if (~isempty(statusHandle))
        % Analysis is complete
        waitbar(1, statusHandle, 'Analysis Complete.');
        close(statusHandle);
    end
    
    % save the parameter file
    parameter_file = fullfile(sesInfo.outputDir, [sesInfo.param_file, '.mat']);
    icatb_save(parameter_file, 'sesInfo');
    
    %% Close parallel mode
    if (strcmpi(parallelMode, 'parallel') && parallelCluster)
        if (~isempty(which('parpool')))
            try
                poolobj = gcp('nocreate');
                delete(poolobj);
            catch
            end
        else
            try
                matlabpool close;
            catch
            end
        end
    end
    
    
    % Use tic and toc instead of cputime
    t_end = toc;
    
    disp(['Time taken to run the analysis is ', num2str(t_end), ' seconds']);
    
    fprintf('\n');
    
    disp(['All the analysis information is stored in the file ', output_LogFile]);
    
    disp('Finished with Analysis');
    fprintf('\n');
    
    diary('off');
    
    if (conserve_disk_space ~= 1)
        analysisComplete = any(stepsToRun == 1) || any(stepsToRun == 7);
    else
        analysisComplete = any(stepsToRun == 1) || any(stepsToRun == 6);
    end
    
    % DO SPM STATS
    if (doSPMStats && analysisComplete)
        % Average runs
        icatb_spm_avg_runs(parameter_file);
        disp('Running one sample t-test on subject component maps ...');
        fprintf('\n');
        
        if spm_stats_avg_runs
            group = (1:sesInfo.numOfSub);
        else
            group = (1:sesInfo.numOfSub*sesInfo.numOfSess);
        end
        
        % Compute one sample t-test on subject maps using SPM
        icatb_spm_stats(parameter_file, 1, 'Group 1', group, [], [], ...
            (1:sesInfo.numComp), spm_stats_write_tal, spm_stats_avg_runs);
        fprintf('\n');
        disp('Done running one sample t-test on subject component maps');
        fprintf('\n');
    end
    % END FOR DOING SPM STATS
    
    if (analysisComplete && (conserve_disk_space == 2))
        cleanupFiles(sesInfo);
    end
    
    
    if open_display_gui
        if (analysisComplete)
            if (strcmpi(modalityType, 'fmri') || strcmpi(modalityType, 'smri'))
                % Open fMRI GUI
                icatb_displayGUI(parameter_file);
            else
                % Open EEG display GUI
                icatb_eeg_displayGUI(parameter_file);
            end
            
        end
    end
    
catch
    
    if exist('sesInfo', 'var')
        % display information to the user
        diary('off');
        if ~strcmpi(analysisType, 'batch')
            % close the status bar
            if (exist('statusHandle', 'var') && ishandle(statusHandle))
                delete(statusHandle);
            end
        end
    end
    
    % Display error message
    icatb_displayErrorMsg;
    
end


function cleanupFiles(sesInfo)
%% Cleanup intermediate analysis files
%

disp('Cleaning up intermediate analysis files ....');
deleteFiles([sesInfo.data_reduction_mat_file, '*.mat'], sesInfo.outputDir, 1);
deleteFiles([sesInfo.back_reconstruction_mat_file, '*.mat'], sesInfo.outputDir, 1);
if (~strcmpi(icatb_get_modality, 'eeg'))
    deleteFiles([sesInfo.calibrate_components_mat_file, '*.mat'], sesInfo.outputDir);
end
disp('Done');
fprintf('\n');

function deleteFiles(fileP, outDir, removeDir)
%% Delete files
%

if (~exist('removeDir', 'var'))
    removeDir = 0;
end

try
    
    p = fileparts(fileP);
    delete(fullfile(outDir, fileP));
    if (~isempty(p) && removeDir)
        rmdir(fullfile(outDir, p));
    end
    
catch
end