function [selectedRegressors, selectedValue] = icatb_returnRegressor(refInfo, num_DataSets, numSubjects, ...
    numSessions, spmMatFlag, subject_string, flagRegressor)
% returns the selected regressor names and values
%%%%%%%%%%%%%%%%%%% DEFINE REGRESSOR GUI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% case 1: datasets = 1, regressors = 1 (no need to select the regressors)
% replicate the regressor over all data-sets
% case 2: datasets = 1, regressors > 1
% replicate the selected regressor over all data-sets
% case 3: datasets > 1, regressors = 1 (no need to select the regressors)
% case 4: datasets > 1, regressors > 1
% depending on the question being asked open the required figure window

% case a: same regressors for all data-sets
% open the list dialog box and select the required regressor and
% replicate over data-sets
% case b: different regressors over sessions
% case 1: number of subjects > 1 & number of sessions = 1
% open the required session regressors in a list dialog box
% else
% open the regressor GUI to pull the session specific
% regressors
% case 2: number of subjects = 1 & number of sessions > 1
% case 3: number of subjects > 1 & number of sessions > 1
% case c: different regressors over subjects & sessions
% open the regressor GUI and specify number of subjects & sessions

%%%%%%%%%%%%%%% END OF DEFINING THE SELECTION OF REGRESSOR %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

icatb_defaults;

global METHOD_ENTERING_REGRESSORS;

if ~strcmpi(METHOD_ENTERING_REGRESSORS, 'batch')
    handle_visiblity = 'on';
else
    handle_visiblity = 'off';
end

if ~exist('flagRegressor', 'var')
    flagRegressor = 'event_avg';
end

% D(1).string = ['For adjusting ICA time-course, only one reference function is used. Please click on the left top listbox to select the data set and the reference functions for that data set will be displayed in the left bottom listbox.'];
% D(size(D,2)+1).string = '';
% D(size(D,2)+1).string = ['Select the reference function in the left bottom listbox and the selected string will appear in the right top listbox. These selections can be made as long as Ok button is pressed.'];
% str = {D.string};

flagRegressor = lower(flagRegressor);

if strcmp(flagRegressor, 'event_avg')
    helpDetails.title = 'Reference Function for Event Avg';
    str = {'For event average of ICA time-course, only one reference function is used. Please select the reference function.'};
    promptStr = ['Reference function for event avg'];
else
    helpDetails.title = 'Adjust ICA time course';
    str = {'For adjusting ICA time-course, only one reference function is included. The parameters other than the selected reference function will be removed.'};
    promptStr = ['Reference function to include in the ICA time course'];
end

% case 1: data sets equal to and the regressor is 1
if num_DataSets  == 1 & size(refInfo.modelIndex, 2) == 1
    % get the selected regressor name

    % selected regressor
    selectedRegressors = deblank(refInfo.selectedRegressors.name);
    
    selectedValue = 1;

elseif num_DataSets  == 1 & size(refInfo.modelIndex, 2) > 1
    % open the list dialog box and select the regressor of interest

    % regressors from sorting GUI
    refNames = str2mat(refInfo.selectedRegressors.name);
    [normalPos] = icatb_getfigPosition('normal');
    % width  and height of the list box
    listsize = [0.75*normalPos(3) 0.9*normalPos(4)];       
        
    helpDetails.str = str;
    
    if ~strcmpi(METHOD_ENTERING_REGRESSORS, 'batch')
        selIndex = icatb_listdlg('PromptString', promptStr, ...
            'SelectionMode','single', 'ListString', refNames, 'listsize', listsize, 'movegui', 'center', ...
            'windowstyle', 'modal', 'help', helpDetails);
    else
        selIndex = 1;
    end
    
    
    if isempty(selIndex)
        disp('Not calculating event average.');
        return;
    end
    
    % selected regressor
    selectedRegressors = deblank(refNames(selIndex, :));
    % selected value
    selectedValue = selIndex;

elseif num_DataSets  > 1 & size(refInfo.modelIndex, 2) == 1
    % regressors from sorting GUI
    selectedRegressors = str2mat(refInfo.selectedRegressors.name);
    % return selected value
    selectedValue = ones(1, size(selectedRegressors, 1));
    
    % check the special case
    if strcmp(lower(spmMatFlag), 'same_sub_same_sess')
        selectedRegressors = repmat(selectedRegressors, num_DataSets, 1);
        % replicate over subjects
        selectedValue = repmat(selectedValue, 1, num_DataSets);
    elseif strcmp(lower(spmMatFlag), 'same_sub_diff_sess')
        % replicate over subjects
        selectedRegressors = repmat(selectedRegressors, numSubjects, 1);
        % replicate over subjects
        selectedValue = repmat(selectedValue, 1, numSubjects);
    end
    
%     % case for session specific regressors
%     if size(selectedRegressors, 1) < numSubjects*numSessions
%         % replicate over subjects
%         selectedRegressors = repmat(selectedRegressors, numSubjects, 1);
%         % replicate over subjects
%         selectedValue = repmat(selectedRegressors, 1, numSubjects);
%     end

elseif num_DataSets  > 1 & size(refInfo.modelIndex, 2) > 1
    % The regressor selection will depend on the answer being chosen to the
    % question being asked in display GUI.

    % check the answer to question in display GUI
    if size(refInfo.SPMFile, 2) == 1 & strcmp(lower(spmMatFlag), 'same_sub_same_sess')
        % this case means one SPM.mat is chosen and the same regressors are
        % applied over all data-sets
        % regressors from sorting GUI
        refNames = str2mat(refInfo.selectedRegressors.name);
        [normalPos] = icatb_getfigPosition('normal');
        % width  and height of the list box
        listsize = [0.75*normalPos(3) 0.9*normalPos(4)];
        %str = {'For event average of ICA time-course, only one reference
        %function is used. Please select the reference function.'};
        helpDetails.str = str;
        
        
        if ~strcmpi(METHOD_ENTERING_REGRESSORS, 'batch')

            selIndex = icatb_listdlg('PromptString', promptStr, ...
                'SelectionMode','single', 'ListString', refNames, 'listsize', listsize, 'movegui', 'center', ...
                'windowstyle', 'modal', 'help', helpDetails);
        else
            selIndex = 1;
        end

        if isempty(selIndex)
            disp('Not calculating event average.');
            return;
        end
        % selected regressor
        selectedRegressors = deblank(refNames(selIndex, :));
        % selected value
        selectedValue = selIndex;
        % replicate over data-sets
        selectedRegressors = repmat(selectedRegressors, num_DataSets, 1);
        % replicate over data-sets
        selectedValue = repmat(selectedValue, 1, num_DataSets);

    elseif size(refInfo.SPMFile, 2) == 1 & strcmp(lower(spmMatFlag), 'same_sub_diff_sess')
        % regressors are different over sessions and one SPM.mat is chosen

        % case b: different regressors over sessions
        % case 1: number of subjects > 1 & number of sessions = 1
        % open the required session regressors in a list dialog box
        % else
        % open the regressor GUI to pull the session specific
        % regressors
        % case 2: number of subjects = 1 & number of sessions > 1
        % case 3: number of subjects > 1 & number of sessions > 1

        if numSubjects > 1 & numSessions == 1
            % regressors from sorting GUI
            refNames = str2mat(refInfo.selectedRegressors.name);
            % number of regressors per session
            numberofRegressors_ses = size(refInfo.modelIndex, 2);
            % pull only the first n regressors
            refNames = refNames(1:numberofRegressors_ses, :);
            [normalPos] = icatb_getfigPosition('normal');
            % width  and height of the list box
            listsize = [0.75*normalPos(3) 0.9*normalPos(4)];
            %str = {'For event average of ICA time-course, only one
            %reference function is used. Please select the reference function.'};
            helpDetails.str = str;
            
            
            if ~strcmpi(METHOD_ENTERING_REGRESSORS, 'batch')

                selIndex = icatb_listdlg('PromptString', promptStr, ...
                    'SelectionMode','single', 'ListString', refNames, 'listsize', listsize, 'movegui', 'center', ...
                    'windowstyle', 'modal', 'help', helpDetails);
            else
                selIndex = 1;
            end

            
            if isempty(selIndex)
                disp('Not calculating event average.');
                return;
            end
            
            % selected regressor
            selectedRegressors = deblank(refNames(selIndex, :));
            % selected value
            selectedValue = selIndex;
            % replicate over data-sets
            selectedRegressors = repmat(selectedRegressors, num_DataSets, 1);
            % replicate over data-sets
            selectedValue = repmat(selectedValue, 1, num_DataSets);

        else
            refNames = str2mat(refInfo.selectedRegressors.name);
            %refNames = refNames(size(refInfo.modelIndex, 2)*numSessions, :);
            % open the session specific regressor GUI
            [refData] = icatb_construct_refData(refNames, numSessions, size(refInfo.modelIndex, 2));
            clear subject_string;
            % construct subject string for sessions
            for ii = 1:numSessions
                subject_string(ii).name = ['Session ', num2str(ii)];
            end
                       
            
            % get the reference data
            refAppData = icatb_selectRef_EventAvg('substring', subject_string, 'numSubjects', 1, ...
                'numSessions', numSessions, 'regressor_data', refData, 'title_figure', promptStr, ...
                'handle_visibility', handle_visiblity);                                            
            
            % selected regressors over sessions
            selectedRegressors = str2mat(refAppData.name);
            % selected value
            selectedValue = [refAppData.value];

            % handling the special case for all subjects and sessions
            if numSubjects > 1 & numSessions > 1
                % replicate over subjects
                selectedRegressors = repmat(selectedRegressors, numSubjects, 1);
                % replicate over subjects
                selectedValue = repmat(selectedValue, 1, numSubjects);
            end
        end

    else
        % add case for all subjects and sessions
        % open figure window for selecting reference functions for calculating
        % event average.
        [refData] = icatb_construct_refData(str2mat(refInfo.selectedRegressors.name), num_DataSets, ...
            size(refInfo.modelIndex, 2));
        refAppData = icatb_selectRef_EventAvg('substring', subject_string, 'numSubjects', ...
            numSubjects, 'numSessions', numSessions, 'regressor_data', refData, 'title_figure', promptStr, ...
            'handle_visibility', handle_visiblity);
        selectedRegressors = str2mat(refAppData.name); % store the sel
        % return selected values
        selectedValue = [refAppData.value];
    end
    % end for checking answer to question in display GUI

end
% end for selecting regressors