function icatb_parRemoveArtifact_cluster(parameter_file, outputDir, selSubjects, selSessions, selComp)
% Remove components from the selected datasets
%
% Inputs:
% 1. parameter_file - Parameter file
% 2. selSubjects - Selected subjects
% 3. selSessions - Selected sessions
% 4. compNum - Component numbers to remove from the data
%


% Defaults
icatb_defaults;
global PARAMETER_INFO_MAT_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
if ~exist('parameter_file', 'var') || isempty(parameter_file)
    parameter_file = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
        filterP, 'title', 'Select a valid parameter file');
end

if isempty(parameter_file)
    error('Parameter file is not selected');
end

if (ischar(parameter_file))
    
    % Results directory
    resultsDir = fileparts(parameter_file);
    
    if isempty(resultsDir)
        resultsDir = pwd;
    end
    
    
    load(parameter_file);
    
else
    sesInfo = parameter_file;
    resultsDir = sesInfo.outputDir;
end

if ~exist('sesInfo', 'var')
    error('Not a valid parameter file. Please select the parameter .MAT file that contains variable sesInfo.');
end

if (~sesInfo.isInitialized)
    error('Please run the analysis in order to remove the components from the data');
end

conserve_disk_space = 0;
if (isfield(sesInfo, 'conserve_disk_space'))
    conserve_disk_space = sesInfo.conserve_disk_space;
end

[modalityType, dataTitle, compSetFields] = icatb_get_modality;


if isfield(sesInfo, 'modality')
    if ~strcmpi(sesInfo.modality, modalityType)
        if strcmpi(sesInfo.modality, 'fmri')
            error('You have selected the fMRI parameter file. Use GIFT toolbox to remove components from the data.');
        elseif strcmpi(sesInfo.modality, 'smri')
            error('You have selected the sMRI parameter file. Use SBM toolbox to remove components from the data.');
        else
            error('You have selected the EEG parameter file. Use EEGIFT toolbox to remove components from the data.');
        end
    end
end


if ~exist('outputDir', 'var') || isempty(outputDir)
    % Select output directory to save images
    outputDir = icatb_selectEntry('typeEntity', 'directory', 'title', 'Select output directory to save images', ...
        'startPath', resultsDir);
end

if isempty(outputDir)
    error('Output directory is not selected');
end

% Get information from the parameter file
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;
numComp = sesInfo.numComp;
inputFiles = sesInfo.inputFiles;
mask_ind = sesInfo.mask_ind;


if (~exist('selSubjects', 'var')) || (isempty(selSubjects))
    % Select subjects
    if numOfSub > 1
        subStr = repmat(struct('name', []), 1, numOfSub);
        for nSub = 1:numOfSub
            subStr(nSub).name = ['Subject ', num2str(nSub)];
        end
        % Select subjects
        [selSubjects] = icatb_listdlg('PromptString', 'Select subjects', 'SelectionMode', ...
            'multiple', 'ListString', char(subStr.name), 'movegui', 'center', ...
            'title_fig', 'Select subjects');
        if isempty(selSubjects)
            error('Subjects need to be selected in order to remove components from the data');
        end
        
    else
        selSubjects = 1;
    end
    % End for selecting subjects
end


if (~exist('selSessions', 'var')) || (isempty(selSessions))
    % Select sessions
    if numOfSess > 1
        sessStr = repmat(struct('name', []), 1, numOfSess);
        for nSess = 1:numOfSess
            sessStr(nSess).name = ['Session ', num2str(nSess)];
        end
        % Select sessions
        [selSessions] = icatb_listdlg('PromptString', 'Select sessions', 'SelectionMode', ...
            'multiple', 'ListString', char(sessStr.name), 'movegui', 'center', ...
            'title_fig', 'Select sessions');
        if isempty(selSessions)
            error('Sessions need to be selected in order to remove components from the data');
        end
    else
        selSessions = 1;
    end
    % End for selecting sessions
    
end

if max(selSubjects) > numOfSub
    error('Error:Num_subjects', ['Check the selSubjects variable as maximum (%s) of selected subjects \nexceeds the number ', ...
        'of subjects (%s)'], num2str(max(selSubjects)), num2str(numOfSub));
end


if max(selSessions) > numOfSess
    error('Error:Num_sessions', ['Check the selSessions variable as maximum (%s) of selected sessions \nexceeds the number ', ...
        'of sessions (%s)'], num2str(max(selSessions)), num2str(numOfSess));
end

if (~exist('selComp', 'var')) || (isempty(selComp))
    % Select components
    compStr = repmat(struct('name', []), 1, numComp);
    for nComp = 1:numComp
        compStr(nComp).name = ['Component ', num2str(nComp)];
    end
    % Select components
    [selComp] = icatb_listdlg('PromptString', 'Select components to remove from the data', 'SelectionMode', ...
        'multiple', 'ListString', char(compStr.name), 'movegui', 'center', ...
        'title_fig', 'Select components');
    if isempty(selComp)
        error('Select components to remove from the data');
    end
    % End for selecting components
end

if max(selComp) > numComp
    error('Error:Num_sessions', ['Check the selComp variable as maximum (%s) of selected components \nexceeds the number ', ...
        'of components (%s)'], num2str(max(selComp)), num2str(numComp));
end

drawnow;

cd(resultsDir);

fprintf('\n');

if (~strcmpi(modalityType, 'smri'))
    disp(['Selected subjects: ', num2str(selSubjects)]);
    disp(['Selected sessions: ', num2str(selSessions)]);
end

% Aggregate spatial maps
try
    load(fullfile(resultsDir, [sesInfo.ica_mat_file, '.mat']), 'icasig');
    icasig = icatb_remove_mean(icasig');
catch
end


numSubjects = length(selSubjects);
numSessions = length(selSessions);

selDataSets = zeros(1, numSubjects*numSessions);

nSet = 0;
for nSub = 1:numSubjects
    for nSess = 1:numSessions
        nSet = nSet + 1;
        selDataSets(nSet) = (selSubjects(nSub) - 1)*numOfSess + selSessions(nSess);
    end
end

back_reconstruction_mat_file = sesInfo.back_reconstruction_mat_file;

disp(['Selected components: ', num2str(selComp)]);


disp(' ');
disp('---------------------------------------------------------------------');
disp('Running component removal tool');
disp('---------------------------------------------------------------------');


parfor nDataSet = 1:length(selDataSets)
    
    nSet = selDataSets(nDataSet);
    
    tmpCompSetFields = compSetFields;
    
    nSub = ceil(nSet/numOfSess);
    nSess = mod(nSet - 1, numOfSess) + 1;
    
    InputFiles = inputFiles;
    
    % Input files
    files = char(InputFiles(nSet).name);
    % Load back-reconstruction MAT file
    subFile = [back_reconstruction_mat_file, num2str(nSet), '.mat'];
    subFile = fullfile(resultsDir, subFile);
    
    loadBackRecon = 0;
    try
        
        if (~conserve_disk_space)
            dd = load(subFile);
            compSet = dd.compSet;
        else
            compSet = getBackReconSet(sesInfo, nSet);
        end
        % Mixing matrix
        A = compSet(1).(tmpCompSetFields{2});
        
        % Force mixing matrix to be timepoints by components
        if size(A, 2) ~= numComp
            A =  A';
        end
        
        loadBackRecon = 1;
        
    catch
        
    end
    
    if (~strcmpi(modalityType, 'smri'))
        % Load data
        msgString = ['Loading Subject ', num2str(nSub), ' Session ', num2str(nSess), ' ...'];
    else
        msgString = 'Loading data ...';
    end
    
    disp(msgString);
    % compute data
    newX = icatb_read_data(files, [], mask_ind);
    % compute the mean
    meanData = mean(newX);
    % remove the mean
    newX = icatb_remove_mean(newX, 0);
    
    if (~loadBackRecon)
        A = (pinv(icasig)*newX)';
    end
    
    % convert to row form (components by volume)
    newX = newX';
    
    S = pinv(A)*newX;
    
    A_R = A; S_R = S;
    A_R(:, selComp) = []; S_R(selComp, :) = [];
    
    
    % compute the new data (contains the removed articfactual)
    disp('Removing the selected components from the data ...');
    newX = newX + A_R*S_R - A*S;
    %clear dewhiteM A S A_R S_R;
    disp('Adding mean back to the data ...');
    % add the mean
    for ii = 1:size(newX, 1)
        newX(ii, :) = newX(ii, :) + meanData(ii);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Step 5: Write new set of images
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%% Create directories %%%%%%
    % Subject directory
    
    subjectDir = icatb_returnFileIndex(nSub);
    
    subjectDir = ['Sub_', subjectDir];
    
    sessDir = num2str(nSess);
    
    if exist(fullfile(outputDir, subjectDir, sessDir)) ~= 7
        mkdir(outputDir, fullfile(subjectDir, sessDir));
    end
    
    sessDir = fullfile(outputDir, subjectDir, sessDir);
    
    msgString = 'Writing out new data-sets ... ';
    
    disp(msgString);
    
    % fMRI, sMRI modality
    if ~strcmpi(modalityType, 'eeg')
        
        % get volume of images
        V = icatb_spm_vol(files);
        % Check for 4d nifti file
        uniqueNames = unique(cellstr(char(V.fname)));
        
        if (length(V) > 1) && (length(uniqueNames) == 1)
            % Get only first file name
            [pathstr, fName, extn] = fileparts(V(1).fname);
            % Add R_ prefix
            newName = ['R_', fName, extn];
            newName = fullfile(sessDir, newName);
            % Read data
            data = icatb_read_vols(V);
            data = reshape(data, [prod(V(1).dim(1:3)), size(data, 4)]);
            % Fill data in regions that contain mask
            data(mask_ind, :) = newX';
            %clear newX;
            descrip = ['Removed component(s) ', num2str(selComp)];
            % Convert back to 4D data
            data = reshape(data, [V(1).dim(1:3), size(data, 2)]);
            % Write Nifti data
            icatb_write_nifti_data(newName, V, data, descrip);
            % clear V data;
        else
            %%%%%%%%% Write out the new images %%%%%%%%%%%%%
            % loop over the number of rows
            for ii = 1:length(V);
                origFileName = V(ii).fname;
                [pathstr, fName, extn] = fileparts(origFileName);
                newName = ['R_', fName, extn];
                V(ii).descrip = ['Removed component(s) ', num2str(selComp)];
                newName = fullfile(sessDir, newName);
                % load the data
                data = icatb_spm_read_vols(V(ii));
                % reshape data
                data = reshape(data, 1, prod(V(ii).dim(1:3)));
                V(ii).fname = newName;
                data(1, mask_ind) = newX(ii, :);
                % reshape the data
                data = reshape(data, V(ii).dim(1:3));
                % write out the data
                icatb_spm_write_vol(V(ii), data);
                %clear data;
            end
            
            %clear V newX;
        end
        
    else
        % EEG modality
        data = icatb_loadData(files);
        size_data = [size(data, 1), size(data, 2), size(data, 3), size(data, 4)];
        data = reshape(data, [prod(size_data(1:3)), size_data(4)]);
        data(mask_ind, :) = newX';
        %clear newX;
        data = reshape(data, size_data);
        
        % Save new data-sets with R_ prefix
        [pathstr, fName, extn] = fileparts(deblank(files(1, :)));
        newName = ['R_', fName, '.mat'];
        newName = fullfile(sessDir, newName);
        icatb_parSave(newName, {data}, {'data'});
        
    end
    % End for checking modalities
    
    msgString = ['The new data-sets are stored in the directory ', sessDir];
    disp(msgString);
    
    fprintf('\n');
    
end
% End loop over number of data-sets

disp('---------------------------------------------------------------------');
disp('Done removing components');
disp('---------------------------------------------------------------------');
disp(' ');

function compSet = getBackReconSet(sesInfo, nDataSet)
%% Get back-reconstructed component of the subject
%

sesInfo.dataSetNo = nDataSet;
[sesInfo, compSet] = icatb_backReconstruct(sesInfo);
