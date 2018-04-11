function [estCompVec, estimateComp, mdlVal, aicVal, sesInfo] = icatb_estimateCompCallback(handleObj, handles, sesInfo )
% Function for doing dimensionality estimation. This function is now
% modified to accept a user specified mask.
%

% input arguments
if ~exist( 'sesInfo', 'var' ) || isempty( sesInfo )
    %  estimate number of components?
    sesInfo = [];
end;
if ( isempty( sesInfo ))
    doGetPrompts = true;
else
    doGetPrompts = false;
end;

% output arguments
if ( nargout > 0 )
    estCompVec = [];
end
if ( nargout > 1 )
    mdlVal = [];
end
if ( nargout > 2 )
    aicVal = [];
end

estimateComp  = 0;

%precisionType = icatb_checkPrecision('single', 0);
precisionType = 'double'; % force double as eigen value solvers can be very sensitive to precision


parallelMode = 'serial';


icatb_defaults;

try
    % get the figure data
    if ( doGetPrompts )
        %sesInfo = get(handles, 'userdata');
        handles_data = get(handles, 'userdata');
        sesInfo = handles_data.sesInfo;
    else
        answerQuestion = 'yes';    % want to estimate number of components?
        answerQuestion2 = 1;       % yes, use all subjects
    end;
    
    output_LogFile = fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, '_results.log']);
    
    % Print output to a file
    diary(output_LogFile);
    
    [sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, 'read');
    outputDir = sesInfo.userInput.pwd;
    dataType = 'real';
    
    icatb_chk_changes_in_mask(sesInfo);
    
    %     maskFile = [];
    %
    %     % mask file
    %     if isfield(sesInfo.userInput, 'maskFile')
    %         maskFile = sesInfo.userInput.maskFile;
    %     end
    
    % check the data type
    if isfield(sesInfo.userInput, 'dataType')
        dataType = sesInfo.userInput.dataType;
    end
    
    if isfield(sesInfo.userInput, 'diffTimePoints')
        diffTimePoints = sesInfo.userInput.diffTimePoints; % store the different time points information
    else
        diffTimePoints = zeros(1, length(sesInfo.userInput.files));
        for ii = 1:length(sesInfo.userInput.files)
            diffTimePoints(ii) = size(sesInfo.userInput.files(ii).name, 1);
        end
        sesInfo.userInput.diffTimePoints = diffTimePoints; % store the different time points information
    end
    
    checkTimePoints = find(diffTimePoints == diffTimePoints(1));
    if length(checkTimePoints) == length(diffTimePoints)
        flagTimePoints = 'same_time_points';
    else
        flagTimePoints = 'different_time_points';
    end
    % check the existence of files
    if ~isfield(sesInfo.userInput, 'files')
        error('Data needs to be selected.');
    end
    
    % Number of subjects
    numOfSub = sesInfo.userInput.numOfSub;
    
    % Number of sessions
    numOfSess = sesInfo.userInput.numOfSess;
    
    %
    if ( doGetPrompts )
        %answerQuestion = get(handleObj, 'value');
        getAnswer = get(handleObj, 'value');
        
        getAnswerString = get(handleObj, 'string');
        
        answerQuestion = getAnswerString(getAnswer, :);
    end;
    
    % If the estimate components is pressed
    if strcmpi(answerQuestion, 'yes') %answerQuestion == 1
        
        % Add question here do you want to select one subject or multiple
        % subjects
        if numOfSub*numOfSess == 1
            
            answerQuestion2 = 1;
            
        else
            
            if ( doGetPrompts )
                questionString1 = 'Do you want to use all subjects to estimate the number of independent components?';
                questionString2 = '';
                questionString3 = ' For every subject the estimated components is calculated and the mean of all the estimated components is taken.';
                questionString4 = '';
                questionString5 = 'Note: When all subjects are loaded it takes lot of computational time.';
                
                questionString = str2mat(questionString1, questionString2, questionString3, questionString4, questionString5);
                
                msgString = cell(1, size(questionString, 1));
                % convert the strings into a cell array
                for i = 1:size(questionString, 1)
                    msgString{i} = questionString(i, :);
                end
                
                [answerQuestion2] = icatb_questionDialog('title', 'Single subject or multiple subjects?', 'textbody', msgString);
                clear questionString questionString1 questionString2 questionString3 questionString4 questionString5 msgString;
            end
            
        end
        
        % Initialise number of components
        estimateComp = 0;
        
        if answerQuestion2 == 1
            
            if (isappdata(0, 'create_mask_gica'))
                % Create mask
                sesInfo = icatb_update_mask(sesInfo);
                rmappdata(0, 'create_mask_gica');
            end
            
            mask = boolean_mask(sesInfo);
            
            num_workers = 4;
            
            try
                parallelMode = sesInfo.userInput.parallel_info.mode;
            catch
            end
            
            try
                num_workers = sesInfo.userInput.parallel_info.num_workers;
            catch
            end
            
            toolboxNames = ver;
            parallelCluster = ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));
            
            if (strcmpi(parallelMode, 'serial'))
                [estimateComp, estCompVec, meanMdl] = serial_computeMDL(sesInfo, handles, mask, precisionType);
            else
                if (parallelCluster)
                    try
                        matlabpool('open', num_workers);
                    catch
                    end
                    [estimateComp, estCompVec, meanMdl] = icatb_parComputeMDL_cluster(sesInfo, mask, precisionType);
                else
                    [estimateComp, estCompVec, meanMdl] = parallel_computeMDL(sesInfo, precisionType);
                end
            end
            
        else
            
            kk = 0;
            questionString = 'Which subject you want to use for estimating components?';
            answerString = cell(numOfSub*numOfSess, 1);
            for ii = 1:numOfSub
                for jj = 1:numOfSess
                    kk = kk + 1;
                    answerString{kk} = ['Subject ', num2str(ii), ' Session ', num2str(jj)];
                end
            end
            
            [answerQuestion3, name_button] = icatb_listdlg('PromptString', questionString, 'SelectionMode','single',...
                'ListString', str2mat(answerString), 'movegui', 'east', 'windowStyle', 'modal');
            
            if(name_button == 0)
                error('Estimating components step is terminated');
                %                 button = msgbox('By default: Selecting the first option', 'Time course is not Selected', 'warn') ;
                %                 waitfor(button);
                %                 answerQuestion3 = 1;
            end
            
            disp(['Loading ', answerString{answerQuestion3}, ' ...']);
            
            if ~isempty(handles)
                helpHandle = helpdlg('System busy. Please wait...', 'Loading  data files');
            end
            % estimate components
            estCompVec = 0;
            drawnow;
            if (isappdata(0, 'create_mask_gica'))
                % Create mask
                sesInfo = icatb_update_mask(sesInfo);
                rmappdata(0, 'create_mask_gica');
            end
            
            mask = boolean_mask(sesInfo);
            
            [estCompVec, mdlVal, aicVal] = icatb_estimate_dimension(sesInfo.userInput.files(answerQuestion3).name, mask, precisionType);
            
            meanMdl = mdlVal;
            
            drawnow;
            %mdl_sub = estCompVec;
            estimateComp = estCompVec;
        end
        try
            if ~isempty(handles)
                if ishandle(helpHandle)
                    delete(helpHandle);
                end
            end
        catch
        end
    end
    
    fprintf('\n');
    if strcmpi(answerQuestion, 'yes') %answerQuestion == 1
        
        % Plot MDL only for one data set
        if answerQuestion2 == 0 || numOfSub*numOfSess == 1
            plotX.y = meanMdl;
            plotX.x = (1:length(meanMdl));
            plotX.title = 'Plot of MDL where minimum is the estimated dimensionality.';
            msgString = ['The estimated independent components is found to be ', num2str(estimateComp), ...
                ' using the MDL criteria.'];
            if ( doGetPrompts )
                helpButton = icatb_dialogBox('title', 'Estimated Components', 'textBody', msgString, 'textType', 'large', 'plotbutton', plotX);
            end
            msgString = {msgString};
        else
            
            meanSub = mean(estCompVec); stdSub = std(estCompVec); [minSub] = min(estCompVec); [maxSub] = max(estCompVec);
            % same time points information
            if strcmp(flagTimePoints, 'same_time_points')
                plotX.y = meanMdl;
                plotX.x = (1:length(meanMdl));
                plotX.title = 'Plot of mean of MDL over subjects.';
            end
            string1 = ['1. The estimated independent components is found to be ', num2str(estimateComp), ...
                ' using the MDL criteria.'];
            string2 = ['2. The mean of mdl over all subjects is found to be ', num2str(meanSub), '.'];
            string3 = ['3. The standard deviation of mdl over all subjects is found to be ', num2str(stdSub), '.'];
            string4 = ['4. The minimum of mdl over all subjects is found to be ', num2str(minSub), '.'];
            string5 = ['5. The maximum of mdl over all subjects is found to be ', num2str(maxSub), '.'];
            string6 = ['6. The median of mdl over all subjects is found to be ', num2str(median(estCompVec)), '.'];
            
            questionString = str2mat(string1, string2, string3, string4, string5, string6);
            msgString = cell(size(questionString, 1), 1);
            % convert the strings into a cell array
            for i = 1:size(questionString, 1)
                msgString{i} = questionString(i, :);
            end
            
            clear questionString;
            
            if ( doGetPrompts )
                if strcmp(flagTimePoints, 'same_time_points')
                    % use the standard dialog box
                    helpButton = icatb_dialogBox('title', 'Estimated Components', 'textBody', msgString, ...
                        'textType', 'large', 'plotbutton', plotX);
                else
                    % use the standard dialog box
                    helpButton = icatb_dialogBox('title', 'Estimated Components', 'textBody', msgString, ...
                        'textType', 'large');
                end
            end;
            
            
        end
        
        for ii = 1:length(msgString)
            fprintf( '%s\n', msgString{ii} );
        end;
        clear msgString;
        
        if (exist('helpButton', 'var'))
            waitfor(helpButton);
        end;
        
    end
    
    if exist('estimateComp', 'var')
        set(handleObj, 'userdata', estimateComp);
    end
    
    sesInfo.userInput.estimated_comps = estCompVec;
    
    diary('off');
catch
    
    diary('off');
    if exist('helpHandle', 'var')
        if ishandle(helpHandle)
            delete(helpHandle);
        end
    end
    
    if ~isempty(handles)
        icatb_errorDialog(lasterr, 'Dimensionality Estimation error');
        icatb_displayErrorMsg;
    else
        icatb_displayErrorMsg;
    end
    
end
% end for estimation

function mask = boolean_mask(sesInfo)

mask = zeros(sesInfo.userInput.HInfo.DIM(1:3));
mask(sesInfo.userInput.mask_ind) = 1;
mask = (mask == 1);


function [estimateComp, estCompVec, meanMdl] = serial_computeMDL(sesInfo, handles, mask, precisionType)


% if (isappdata(0, 'create_mask_gica'))
%     % Create mask
%     sesInfo = icatb_update_mask(sesInfo);
%     rmappdata(0, 'create_mask_gica');
% end
%
% mask = boolean_mask(sesInfo);

estimateComp = 0;

if ~isempty(handles)
    helpHandle = helpdlg('System busy. Please wait...', 'Loading  data files');
end

% Load functional data files here and mask out voxels
% and then call estimate dimension to return the number of components
estCompVec = zeros( length(sesInfo.userInput.files), 1 );
tpMin = min(sesInfo.userInput.diffTimePoints) - 1;
meanMdl = zeros(1, tpMin);
drawnow;
for numFiles = 1:length(sesInfo.userInput.files)
    fprintf('\n');
    disp('..............................................');
    disp('');
    helpHandleName = ['Loading data set ', num2str(numFiles)];
    disp(helpHandleName);
    if ~isempty(handles)
        if ishandle(helpHandle)
            set(helpHandle, 'name', helpHandleName);
        end
    end
    drawnow;
    
    [ estCompVec(numFiles), mdlVal] = icatb_estimate_dimension(sesInfo.userInput.files(numFiles).name, mask, precisionType);
    %     if strcmp(flagTimePoints, 'same_time_points')
    %
    %         if length(sesInfo.userInput.files) > 1
    %             if numFiles == 1
    %                 meanMdl = zeros(size(mdlVal));
    %             end
    %             meanMdl = meanMdl + mdlVal;
    %         end
    %
    %     end
    
    dd = min([tpMin, length(mdlVal)]);
    
    meanMdl(1:dd) = meanMdl(1:dd) + mdlVal(1:dd);
    
    estimateComp = estCompVec(numFiles) + estimateComp;
    
    drawnow;
    
end

meanMdl = meanMdl/length(sesInfo.userInput.files);

estimateComp = round(estimateComp/length(sesInfo.userInput.files));


if ~isempty(handles)
    if ishandle(helpHandle)
        delete(helpHandle);
    end
end


function [estimateComp, estCompVec, meanMdl] = parallel_computeMDL(sesInfo, precisionType)
%% Compute MDL in parallel mode (sessions)
%

cValue = clock;
% remove a minute
cValue(end-1) = max([0, cValue(end-1) - 1]);
currentTime = datenum(cValue);

%% Initialize varaibles
totalSubjectsIn = length(sesInfo.userInput.files);
dataSetsToRun = (1:totalSubjectsIn);
num_workers = 4;
try
    num_workers = sesInfo.userInput.parallel_info.num_workers;
catch
end

if (num_workers > totalSubjectsIn)
    num_workers = totalSubjectsIn;
end

disp('Computing MDL for all subjects and sessions in different MATLAB sessions. Please wait ....');

drawnow;

paramFile = sesInfo.userInput.param_file;
icatb_save(paramFile, 'sesInfo');
giftPath = fileparts(which('gift.m'));
dummyScriptPath = fullfile(giftPath, 'icatb_parallel_files', 'icatb_dummyScript.m');
increments = ceil(totalSubjectsIn/num_workers);
mdlFiles = cellstr(strcat(fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, '_mdl_inc_']), icatb_numberToString((1:num_workers)'), '.mat'));

%% Run MDL in separate sessions
eW = 0;
for nF = 1:num_workers
    sW = eW + 1;
    eW = eW + increments;
    eW = min([totalSubjectsIn, eW]);
    tmpDataSetsToRun = dataSetsToRun(sW:eW);
    % Run separate matlab sessions in background mode
    % (Dummyscript i.e., no text required to run in linux OS)
    commandStr = ['!matlab -nodesktop -nosplash -r "addpath(genpath(''', giftPath, '''));icatb_parComputeMDL(''', paramFile, ''',[', ...
        num2str(tmpDataSetsToRun), '],''', precisionType, ''',''', mdlFiles{nF}, ''');exit" < "', dummyScriptPath, '" &'];
    eval(commandStr);
end

icatb_waitForTaskCompletion(mdlFiles, currentTime);

%% Gather output
tpMin = min(sesInfo.userInput.diffTimePoints) - 1;
meanMdl = zeros(1, tpMin);
estCompVec = [];
for nOut = 1:length(mdlFiles)
    p=load(mdlFiles{nOut});
    meanMdl = meanMdl + p.meanMdl;
    estCompVec = [estCompVec, p.estCompVec];
    delete(mdlFiles{nOut});%cleanup
end

estimateComp = round(mean(estCompVec));
meanMdl = meanMdl/totalSubjectsIn;