function icatb_showTimeCourses(varargin)
% Multiple plots incorporated in this function. Figure contains time courses with or
% without model, with or without mean, with or without SEM (standard error
% of mean).
%
% Input:
% 1. timeCourseStruct
% 2. modelTimeCourseStruct
% 3. xAxis
% 4. plotColors
% 5. timeCourseNames to be plotted as figure title
% 6. meanDisplay ('yes' or no')
% 7. SEMDisplay ('yes' or 'no')
% 8. figure title
% Output:

if nargin == 0
    error('Arguments are needed');
end

if mod(nargin, 2) ~= 0
    error('Arguments must be in pairs');
end

% Initialise vars or defaults
meanDisplay = 'no';
calculateSEM = 'no';
timeCourseColor = 'c';
modelTimeCourseColor = 'y:';
SEM_Color = 'g:';
countTimePoints = [];
modelTimeCourseStruct = [];
flagXAxis = 'yes';
% rows contain subjects for session timecourses
% columns contain subject's session timecourses
% at the end of cols session mean can be taken
% at the end of rows subject mean can be taken

% Loop over number of args
for ii = 1:2:nargin
    % get the time course strcuture
    % time course containing subjects and sessions
    if strcmpi(varargin{ii}, 'tc')
        timeCourseStruct = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'model')
        modelTimeCourseStruct = varargin{ii + 1};
        % mean display ('no' or 'yes')
    elseif strcmpi(varargin{ii}, 'meandisplay')
        meanDisplay = lower(varargin{ii + 1});
        % x-axis structure
    elseif strcmpi(varargin{ii}, 'x_axis');
        x_Axis = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'titlefig');
        % title for the figure
        titleFig = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'sem');
        % title for the figure
        calculateSEM = lower(varargin{ii + 1});
    elseif strcmpi(varargin{ii}, 'difftimepoints');
        countTimePoints = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'timecourse_color')
        timeCourseColor = lower(varargin{ii + 1});
    elseif strcmpi(varargin{ii}, 'model_timecourse_color')
        modelTimeCourseColor = lower(varargin{ii + 1});
    elseif strcmpi(varargin{ii}, 'sem_color')
        SEM_Color = lower(varargin{ii + 1});
    elseif strcmpi(varargin{ii}, 'title_axes')
        tag = varargin{ii + 1};
    end
end

% number of subjects and sessions
numSubjects = length(timeCourseStruct.sub);
numSessions = length(timeCourseStruct.sub(1).sess);

% mean is not calculated if the number of subjects*sessions is 1
if numSessions*numSubjects == 1
    meanDisplay = 'no';
end

if ~exist('x_Axis', 'var')
    flagXAxis = 'no';
end

% check the length of the model time course and ICA time course
if exist('modelTimeCourseStruct', 'var')
    if ~isempty(modelTimeCourseStruct)
        if ~isfield(modelTimeCourseStruct, 'sub')
            error('modelTimeCourseStruct should contains sub field');
        end
        if ~isfield(modelTimeCourseStruct.sub, 'sess')
            error('modelTimeCourseStruct.sub should contains sess field');
        end
        if length(modelTimeCourseStruct.sub) ~= numSubjects
            error('modelTimeCourseStruct subjects doesn''t match the number of subjects');
        end
        if length(modelTimeCourseStruct.sub(1).sess) ~= numSessions
            error('modelTimeCourseStruct sessions doesn''t match the number of sessions');
        end
    end
end

% this function has to handle the different time points
% this function should plot should have the ability to plot model
if isempty(countTimePoints)
    %%%%%%% Checking time points and model time course %%%%%
    % check if all have the same time points or not
    countDataSets = 0;
    countTimePoints = zeros(1, numSubjects*numSessions);
    for numSub = 1:numSubjects
        for numSess = 1: numSessions
            countDataSets = countDataSets + 1;
            % time course
            countTimePoints(countDataSets) = length(timeCourseStruct.sub(numSub).sess(numSess).tc);
            if strcmpi(flagXAxis, 'no')
                x_Axis.sub(numSub).sess(numSess).tc = (1:countTimePoints(countDataSets))';
            end
            if prod(size(timeCourseStruct.sub(numSub).sess(numSess).tc)) ~= countTimePoints(countDataSets)
                error(['Timecourse for subject ', num2str(numSub), ' session ', num2str(numSess), ...
                    ' should be a vector']);
            end
            % end for checking time course
            if ~isempty(modelTimeCourseStruct)
                if length(modelTimeCourseStruct.sub(numSub).sess(numSess).tc) ~= countTimePoints(countDataSets)
                    error(['Model time course for subject ', num2str(numSub), ' session ', num2str(numSess), ...
                        ' should be same as the ICA time course']);
                end
                if size(modelTimeCourseStruct.sub(numSub).sess(numSess).tc, 1) ~= countTimePoints(countDataSets)
                    currentTc = modelTimeCourseStruct.sub(numSub).sess(numSess).tc;
                    modelTimeCourseStruct.sub(numSub).sess(numSess).tc = currentTc';
                end
            end
            % end for checking model time course
        end
    end
    % end for loop over data sets

end

checkTimePoints = find(countTimePoints ~= countTimePoints(1) );

nTimePoints =  countTimePoints(1); % calculate from the first time course

clear countTimePoints;

% don't calculate mean if time points are different
if ~isempty(checkTimePoints)
    meanDisplay = 'no';
end

% % don't calculate sem if mean display is off
if strcmpi(meanDisplay, 'no')
    calculateSEM = 'no';
end

% Initialise variables
meanAll_sub_sess = []; meanSessionTc = []; meanSubjectTc = [];

% mean display is not off
% compute mean over rows, cols and mean of rows and cols
if ~strcmpi(meanDisplay, 'no')
    helpHandle = helpdlg('Calculating stats. Please Wait...');
    % Initialise SEM_sessions
    SEM_sessions = zeros(nTimePoints, numSubjects);
    % Initialise SEM_subjects
    SEM_subjects = zeros(nTimePoints, numSessions);
    % compute mean over subject's sessions
    meanSessionTc = zeros(nTimePoints, numSubjects);
    for numSub = 1:numSubjects
        sessionTc = zeros(nTimePoints, numSessions);
        for numSess = 1:numSessions
            sessionTc(:, numSess) = timeCourseStruct.sub(numSub).sess(numSess).tc;
        end
        meanSessionTc(:, numSub) = mean(sessionTc, 2);

        if ~strcmpi(calculateSEM, 'no') & numSessions ~= 1
            tempSess = std(sessionTc, 0, 2);
            SEM_sessions(:, numSub) = tempSess./sqrt(numSessions - 1);
        end
        clear sessionTc;
        clear tempSess;
    end

    meanSubjectTc = zeros(nTimePoints, numSessions);
    % compute mean over subjects
    for numSess = 1:numSessions
        subjectTc = zeros(nTimePoints, numSubjects);
        for numSub = 1:numSubjects
            subjectTc(:, numSub) = timeCourseStruct.sub(numSub).sess(numSess).tc;
        end
        meanSubjectTc(:, numSess) = mean(subjectTc, 2);

        if ~strcmpi(calculateSEM, 'no') & numSubjects ~= 1
            tempSub = std(subjectTc, 0, 2);
            SEM_subjects(:, numSess) = tempSub./sqrt(numSubjects - 1);
        end

        clear subjectTc;
        clear tempSub;
    end
    % mean over all subjects and sessions
    meanAll_sub_sess = mean(meanSubjectTc, 2);

    if ~strcmpi(calculateSEM, 'no') & numSubjects ~= 1 & numSessions ~= 1
        meanAll_SEM = zeros(nTimePoints, 1);
        %meanAll_SEM = mean(SEM_subjects, 2);
        allTimeCourses = zeros(nTimePoints, numSubjects*numSessions);
        count = 0;
        for allSub = 1:numSubjects
            for allSess = 1:numSessions
                count = count + 1;
                allTimeCourses(1:nTimePoints, count) = timeCourseStruct.sub(allSub).sess(allSess).tc;
            end
        end
        meanAll_SEM = std(allTimeCourses, 0, 2)./ sqrt((numSubjects*numSessions) - 1);
        clear allTimeCourses;
    end
    try
        close(helpHandle);
    catch
    end
end
% end for group stats

helpHandle = helpdlg('Displaying Time courses. Please Wait...', 'Displaying Time courses');
% load defaults file
icatb_defaults;
global FONT_COLOR;
global BG_COLOR;
global BG2_COLOR;


%helpHandle = helpdlg('Plotting time courses. Please Wait...');

if ~exist('titleFig', 'var')
    titleFig = 'Figure showing plots of multiple subjects and sessions';
    set(graphicsHandle, 'tag', titleFig);
end

% Set the figure properties here
graphicsHandle = icatb_getGraphics(titleFig, 'Graphics', titleFig);


% GUI Options Menu
OptionsMenu = uimenu('parent', graphicsHandle, 'label', 'Options');
extract_tc_Menu = uimenu(OptionsMenu, 'label', 'Extract ICA Timecourses');

keywd = 'Subjects and Sessions';

% axes dimensions
axesWidth = 0.4; axesHeight = 0.4;

% X Offset and YOffset
xOffset = 0.05; yOffset = 0.08;

% Initial Axes Position
axesYOrigin =  1 - yOffset - axesHeight;

temp = axesYOrigin;

% Plot in case of one session for a subject
numFiguresPerCol = 2;

sliderStep = [0.05 0.2];

% 2 cases to deal with
% first case contains either of numSubjects and numSessions == 1
% in the first case condition is also included to incorporate numSubjects
% == 1 and numSessions == 1
%
if (numSessions*numSubjects == 1) ||  (numSessions == 1) || (numSubjects == 1)
    if numSessions*numSubjects == 1
        keywd = 'none';
        % Number of subjects and mean over subjects
        numdataSets = numSubjects + 1;
    end
    if numSessions == 1
        keywd = 'subjects';
        % Number of subjects and mean over subjects
        numdataSets = numSubjects + 1;
    end
    if numSubjects == 1
        keywd = 'sessions';
        numdataSets = numSessions + 1;
    end

    % number of rows
    nrows = ceil(numdataSets / numFiguresPerCol);

    yPositions = zeros(1, nrows);
    % determine the y axis positions if only vertical scroller is used
    for nRow = 1:nrows
        yPositions(nRow) = temp;
        temp = temp - axesHeight - yOffset;
    end

else
    numdataSets = (numSessions + 1)*(numSubjects + 1);
end

if strcmpi(meanDisplay, 'no')
    numdataSets = numSessions*numSubjects;
end

if ~exist('tag', 'var')
    tag = cell(numdataSets, 1);
    if ~strcmpi(meanDisplay, 'no')
        if strcmpi(keywd, 'subjects')
            % Loop over number of datasets
            for nSub = 1:numSubjects + 1
                if nSub ~= numSubjects + 1
                    tag{nSub} = ['Subject ', num2str(nSub)];
                else
                    if ~strcmpi(calculateSEM, 'no')
                        tag{nSub} = 'Mean over subjects + SEM';
                    else
                        tag{nSub} = 'Mean over subjects';
                    end
                end
            end
            % End loop over subjects
        elseif strcmpi(keywd, 'sessions')
            % Loop over number of datasets
            for nSess = 1:numSessions + 1
                if nSess ~= numSessions + 1
                    tag{nSess} = ['Session ', num2str(nSess)];
                else
                    if ~strcmpi(calculateSEM, 'no')
                        tag{nSess} = 'Mean over sessions + SEM';
                    else
                        tag{nSess} = 'Mean over sessions';
                    end
                end
            end
            % End loop over data-sets
        else

            countTag = 0;
            % Loop over subjects + 1
            for nSub = 1:numSubjects + 1
                % Loop over sessions + 1
                for nSess = 1:numSessions + 1
                    countTag = countTag + 1;
                    if  (nSub ~= numSubjects + 1) && (nSess ~= numSessions + 1)
                        tag{countTag} = ['Subject ', num2str(nSub), ' Session ', num2str(nSess)];
                    elseif (nSub ~= numSubjects + 1) && (nSess == numSessions + 1)
                        if ~strcmpi(calculateSEM, 'no')
                            tag{countTag} = ['Subject ', num2str(nSub), ' Session Mean + SEM'];
                        else
                            tag{countTag} = ['Subject ', num2str(nSub), ' Session Mean'];
                        end
                    elseif (nSub == numSubjects + 1) && (nSess ~= numSessions + 1)
                        if ~strcmpi(calculateSEM, 'no')
                            tag{countTag} = ['Mean over Subjects of Session ', num2str(nSess), ' + SEM'];
                        else
                            tag{countTag} = ['Mean over Subjects of Session ', num2str(nSess)];
                        end
                    else
                        if ~strcmpi(calculateSEM, 'no')
                            tag{countTag} = 'Mean over all subjects and sessions + SEM';
                        else
                            tag{countTag} = 'Mean over all subjects and sessions';
                        end
                    end
                end
                % End loop over sessions + 1
            end
            % End loop over subjects + 1

        end
        % End for checking

    else

        if strcmpi(keywd, 'subjects')
            % Loop over number of datasets
            for nSub = 1:numSubjects
                tag{nSub} = ['Subject ', num2str(nSub)];
            end
            % End loop over subjects
        elseif strcmpi(keywd, 'sessions')
            % Loop over number of datasets
            for nSess = 1:numSessions
                tag{nSess} = ['Session ', num2str(nSess)];
            end
            % End loop over data-sets
        else

            countTag = 0;
            % Loop over subjects
            for nSub = 1:numSubjects
                % Loop over sessions
                for nSess = 1:numSessions
                    countTag = countTag + 1;
                    tag{countTag} = ['Subject ', num2str(nSub), ' Session ', num2str(nSess)];
                end
                % End loop over sessions
            end
            % End loop over subjects

        end
        % End for checking

    end
end

if length(tag) ~= numdataSets
    error('Error:Title_Fig', 'Number of titles (%s) doesn''t match no. of data-sets (%s)', num2str(length(tag)), num2str(numdataSets));
end


% Initialise initial positions
initialPosition = zeros(numdataSets, 4);
axesHandles = zeros(numdataSets, 1);

% Plots subject's time courses in 2 columns
if numSessions == 1 || numSubjects == 1
    % loop over figures
    for nDataSets = 1:numdataSets
        % Getting the position in x of the axes
        if mod(nDataSets, numFiguresPerCol) ~= 0
            % Axes X origin
            axesXOrigin = xOffset;
        else
            % Axes X origin
            axesXOrigin = axesWidth + 2*xOffset;
        end
        % get the row number
        rowNumber = ceil(nDataSets/numFiguresPerCol);
        % get the related y position
        axesYOrigin = yPositions(rowNumber);
        % Plot all the individual subjects here
        if ((nDataSets ~= numdataSets) || strcmpi(meanDisplay, 'no'))
            % set axes with tag and position
            axesH = axes('Parent', graphicsHandle, 'units', 'normalized', 'position', ...
                [axesXOrigin, axesYOrigin, axesWidth, axesHeight]);

            % Set tag to the axes corresponding to sessions
            if strcmpi(keywd, 'sessions')
                % store time courses in extract data
                extractData.timecourses(nDataSets).tc = timeCourseStruct.sub(1).sess(nDataSets).tc;
                extractData.xAxis(nDataSets).tc = x_Axis.sub(1).sess(nDataSets).tc;
                % Plot ICA time course
                plot(extractData.xAxis(nDataSets).tc, extractData.timecourses(nDataSets).tc, ...
                    timeCourseColor, 'parent', axesH);
                hold on;
                %tag{nDataSets} = ['Session ', num2str(nDataSets)];
                title(tag{nDataSets}, 'parent', axesH);
                % if model time course is present
                if ~isempty(modelTimeCourseStruct)
                    nRegress = size(modelTimeCourseStruct.sub(1).sess(nDataSets).tc, 2); % number of regressors
                    for num_regress = 1:nRegress
                        hold on;
                        plot(x_Axis.sub(1).sess(nDataSets).tc, modelTimeCourseStruct.sub(1).sess(nDataSets).tc(:, num_regress), ...
                            modelTimeCourseColor, 'parent', axesH);
                    end
                    hold off;
                end
            elseif strcmpi(keywd, 'subjects')
                % store time courses in extract data
                extractData.timecourses(nDataSets).tc = timeCourseStruct.sub(nDataSets).sess(1).tc;
                extractData.xAxis(nDataSets).tc = x_Axis.sub(nDataSets).sess(1).tc;
                % Plot ICA time course
                plot(extractData.xAxis(nDataSets).tc, extractData.timecourses(nDataSets).tc, ...
                    timeCourseColor, 'parent', axesH);
                %tag{nDataSets} = ['Subject ', num2str(nDataSets)];
                title(tag{nDataSets}, 'parent', axesH);
                % if model time course is present
                if ~isempty(modelTimeCourseStruct)
                    nRegress = size(modelTimeCourseStruct.sub(nDataSets).sess(1).tc, 2); % number of regressors
                    for num_regress = 1:nRegress
                        hold on;
                        plot(x_Axis.sub(nDataSets).sess(1).tc, modelTimeCourseStruct.sub(nDataSets).sess(1).tc(:, num_regress), ...
                            modelTimeCourseColor, 'parent', axesH);
                    end
                    hold off;
                end
            else
                % store time courses in extract data
                extractData.timecourses(nDataSets).tc = timeCourseStruct.sub(1).sess(1).tc;
                extractData.xAxis(nDataSets).tc = x_Axis.sub(1).sess(1).tc;
                % Plot ICA time course
                plot(extractData.xAxis(nDataSets).tc, extractData.timecourses(nDataSets).tc, ...
                    timeCourseColor, 'parent', axesH);
                %tag{nDataSets} = ['Subject ', num2str(nDataSets), ' Session ', num2str(nDataSets)];
                title(tag{nDataSets}, 'parent', axesH);
                % if model time course is present
                if ~isempty(modelTimeCourseStruct)
                    nRegress = size(modelTimeCourseStruct.sub(1).sess(nDataSets).tc, 2); % number of regressors
                    for num_regress = 1:nRegress
                        hold on;
                        plot(x_Axis.sub(1).sess(nDataSets).tc, modelTimeCourseStruct.sub(1).sess(nDataSets).tc(:, num_regress), ...
                            modelTimeCourseColor, 'parent', axesH);
                    end
                    hold off;
                end
            end
            axis(axesH, 'tight');
            set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
            set(axesH, 'tag', tag{nDataSets});
            % Store all the initial positions
            initialPosition(nDataSets, :) = [axesXOrigin, axesYOrigin, axesWidth, axesHeight]; %get(axesH, 'position');
            axesHandles(nDataSets) = axesH;
            lastOrigin = axesXOrigin;
            lastYOrigin = axesYOrigin;
        else
            % Plot the mean if the display is on
            if ~strcmpi(meanDisplay, 'no')
                % Plot mean timecourse
                % set axes with tag and position
                meanAxesH = axes('Parent', graphicsHandle, 'units', 'normalized', 'position', ...
                    [axesXOrigin, axesYOrigin, axesWidth, axesHeight]);
                % Set tag to the axes corresponding to sessions
                if strcmpi(keywd, 'sessions')
                    %tag{nDataSets} = ['Mean over sessions'];
                    % store time courses in extract data
                    extractData.timecourses(nDataSets).tc = meanSessionTc;
                    extractData.xAxis(nDataSets).tc = x_Axis.sub(1).sess(1).tc;
                    % Plot ICA time course
                    plot(extractData.xAxis(nDataSets).tc, extractData.timecourses(nDataSets).tc, timeCourseColor, 'parent', meanAxesH);
                    % plot if Calculate Standard Error Mean exists
                    if ~strcmpi(calculateSEM, 'no') & numSessions ~= 1
                        hold on;
                        % Calculate Standard Error mean
                        SEM_meanSessionTc_up = SEM_sessions + meanSessionTc;
                        plot(x_Axis.sub(1).sess(1).tc, SEM_meanSessionTc_up, SEM_Color, 'parent', meanAxesH);
                        hold on;
                        SEM_meanSessionTc_down = meanSessionTc - SEM_sessions;
                        plot(x_Axis.sub(1).sess(1).tc, SEM_meanSessionTc_down, SEM_Color, 'parent', meanAxesH);
                        hold off;
                        clear SEM_meanSessionTc_down; clear SEM_meanSessionTc_up;
                        %tag{nDataSets} = ['Mean over sessions + SEM'];
                    end
                    title(tag{nDataSets}, 'parent', meanAxesH);
                end
                % end for sessions

                % Set tag to the axes corresponding to subjects
                if strcmpi(keywd, 'subjects')
                    %tag{nDataSets} = ['Mean over subjects'];
                    % store time courses in extract data
                    extractData.timecourses(nDataSets).tc = meanSubjectTc;
                    extractData.xAxis(nDataSets).tc = x_Axis.sub(1).sess(1).tc;
                    % Plot ICA time course
                    plot(extractData.xAxis(nDataSets).tc, extractData.timecourses(nDataSets).tc, timeCourseColor, 'parent', meanAxesH);
                    % plot if Calculate Standard Error Mean exists
                    if ~strcmpi(calculateSEM, 'no') & numSubjects ~= 1
                        hold on;
                        % Calculate Standard Error mean
                        SEM_meanSubjectTc_up = SEM_subjects + meanSubjectTc;
                        plot(x_Axis.sub(1).sess(1).tc, SEM_meanSubjectTc_up, SEM_Color, 'parent', meanAxesH);
                        hold on;
                        SEM_meanSubjectTc_down = meanSubjectTc - SEM_subjects;
                        plot(x_Axis.sub(1).sess(1).tc, SEM_meanSubjectTc_down, SEM_Color, 'parent', meanAxesH);
                        hold off;
                        %tag{nDataSets} = ['Mean over subjects + SEM'];
                    end
                    title(tag{nDataSets}, 'parent', meanAxesH);
                end
                % end for subjects
                axis(meanAxesH, 'tight');
                set(meanAxesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
                set(meanAxesH, 'tag', tag{nDataSets});
                % Store all the initial positions
                initialPosition(nDataSets, :) = [axesXOrigin, axesYOrigin, axesWidth, axesHeight]; %get(meanAxesH, 'position');
                axesHandles(nDataSets) = meanAxesH;
                lastOrigin = axesXOrigin;
                lastYOrigin = axesYOrigin;
            end
            % end for mean display
        end
        % end for checking data sets number
    end
    % end for loop over data sets
else
    % loop over all data sets
    % Initialise variables
    count = 0;
    % Loop over subjects
    for numSub = 1 : numSubjects
        % Axes X origin
        axesXOrigin = 0.05;
        % Loop over sessions
        for numSess = 1 : numSessions
            count = count + 1;
            % if the number doesn''t match the mean for the sessions
            % number
            % set axes with tag and position
            axesH = axes('Parent', graphicsHandle, 'units', 'normalized', 'position', ...
                [axesXOrigin, axesYOrigin, axesWidth, axesHeight]);
            %tag{count} = ['Subject ', num2str(numSub), ' Session ', num2str(numSess)];
            % store time courses in extract data
            extractData.timecourses(count).tc = timeCourseStruct.sub(numSub).sess(numSess).tc;
            extractData.xAxis(count).tc = x_Axis.sub(numSub).sess(numSess).tc;
            % Plot ICA time course
            plot(extractData.xAxis(count).tc, extractData.timecourses(count).tc, ...
                timeCourseColor, 'parent', axesH);
            title(tag{count}, 'parent', axesH);
            % if model time course is present
            if ~isempty(modelTimeCourseStruct)
                nRegress = size(modelTimeCourseStruct.sub(numSub).sess(numSess).tc, 2); % number of regressors
                for num_regress = 1:nRegress
                    hold on;
                    plot(x_Axis.sub(numSub).sess(numSess).tc, modelTimeCourseStruct.sub(numSub).sess(numSess).tc(:, num_regress), ...
                        modelTimeCourseColor, 'parent', axesH);
                end
                hold off;
            end
            axis(axesH, 'tight');
            set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
            set(axesH, 'tag', tag{count});
            % Store positions of axes
            initialPosition(count, :) = [axesXOrigin, axesYOrigin, axesWidth, axesHeight]; %get(axesH, 'position');
            axesHandles(count) = axesH;
            lastOrigin = axesXOrigin;
            lastYOrigin = axesYOrigin;
            axesXOrigin = axesXOrigin + axesWidth + xOffset;

            if ~strcmpi(meanDisplay, 'no') & numSess == numSessions
                count = count + 1;
                % plot mean over sessions
                % get the mean of all the time courses
                % set axes with tag and position
                meanAxesH = axes('Parent', graphicsHandle, 'units', 'normalized', 'position', ...
                    [axesXOrigin, axesYOrigin, axesWidth, axesHeight]);
                %tag{count} = ['Subject ', num2str(numSub), ' Session Mean'];
                % store time courses in extract data
                extractData.timecourses(count).tc = meanSessionTc(:, numSub);
                extractData.xAxis(count).tc = x_Axis.sub(1).sess(1).tc;
                % Plot ICA time course
                plot(extractData.xAxis(count).tc, extractData.timecourses(count).tc, timeCourseColor, 'parent', meanAxesH);
                % plot if Calculate Standard Error Mean exists
                if ~strcmpi(calculateSEM, 'no')
                    hold on;
                    % Calculate Standard Error mean
                    SEM_meanSessionTc_up = SEM_sessions(:, numSub) + meanSessionTc(:, numSub);
                    plot(x_Axis.sub(1).sess(1).tc, SEM_meanSessionTc_up, SEM_Color, 'parent', meanAxesH);
                    hold on;
                    SEM_meanSessionTc_down = meanSessionTc(:, numSub) - SEM_sessions(:, numSub);
                    plot(x_Axis.sub(1).sess(1).tc, SEM_meanSessionTc_down, SEM_Color, 'parent', meanAxesH);
                    hold off;
                    %tag{count} = ['Subject ', num2str(numSub), ' Session Mean + SEM'];
                    clear SEM_meanSessionTc_down; clear SEM_meanSessionTc_up;
                end
                title(tag{count}, 'parent', meanAxesH);
                axis(meanAxesH, 'tight');
                set(meanAxesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
                set(meanAxesH, 'tag', tag{count});
                % Store positions of axes
                initialPosition(count, :) = [axesXOrigin, axesYOrigin, axesWidth, axesHeight]; %get(meanAxesH, 'position');
                axesHandles(count) = meanAxesH;

                lastOrigin = axesXOrigin;
                lastYOrigin = axesYOrigin;
            end
            if numSess == numSessions
                axesYOrigin = axesYOrigin - axesHeight - yOffset;
            end
        end
        % end for loop over sessions
    end
    % end for loop over subjects

    % Axes X origin
    axesXOrigin = initialPosition(1, 1);
    % when mean display is not off
    if ~strcmpi(meanDisplay, 'no')
        % Compute mean over subjects
        for numSess = 1 : numSessions
            % Update the count here
            count = count + 1;
            % if the number is not equal to the mean over all subjects and sessions
            % plot mean over subjects
            % set axes with tag and position
            meanAxesH = axes('Parent', graphicsHandle, 'units', 'normalized', 'position', ...
                [axesXOrigin, axesYOrigin, axesWidth, axesHeight]);
            % store time courses in extract data
            extractData.timecourses(count).tc = meanSubjectTc(:, numSess);
            extractData.xAxis(count).tc = x_Axis.sub(1).sess(1).tc;
            % Plot ICA time course
            plot(extractData.xAxis(count).tc, extractData.timecourses(count).tc, timeCourseColor, 'parent', meanAxesH);
            %tag{count} = ['Mean over Subjects of Session ', num2str(numSess)];
            % plot if Calculate Standard Error Mean exists
            if ~strcmpi(calculateSEM, 'no')
                hold on;
                % Calculate Standard Error mean
                SEM_meanSubjectTc_up = SEM_subjects(:, numSess) + meanSubjectTc(:, numSess);
                plot(x_Axis.sub(1).sess(1).tc, SEM_meanSubjectTc_up, SEM_Color, 'parent', meanAxesH);
                hold on;
                SEM_meanSubjectTc_down = meanSubjectTc(:, numSess) - SEM_subjects(:, numSess);
                plot(x_Axis.sub(1).sess(1).tc, SEM_meanSubjectTc_down, SEM_Color, 'parent', meanAxesH);
                hold off;
                clear SEM_meanSubjectTc_up;
                clear SEM_meanSubjectTc_down;
                %tag{count} = ['Mean over Subjects of Session ', num2str(numSess), ' + SEM'];
            end
            title(tag{count}, 'parent', meanAxesH);
            axis(meanAxesH, 'tight');
            set(meanAxesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
            set(meanAxesH, 'tag', tag{count});
            % Store all positions
            initialPosition(count, :) = [axesXOrigin, axesYOrigin, axesWidth, axesHeight]; %get(meanAxesH, 'position');
            axesHandles(count) = meanAxesH;
            lastOrigin = axesXOrigin;
            lastYOrigin = axesYOrigin;
            axesXOrigin = axesXOrigin + axesWidth + xOffset; % update axes x position
            if numSess == numSessions
                count = count + 1;
                % set axes with tag and position
                meanAxesH = axes('Parent', graphicsHandle, 'units', 'normalized', 'position', ...
                    [axesXOrigin, axesYOrigin, axesWidth, axesHeight]);
                % store time courses in extract data
                extractData.timecourses(count).tc = meanAll_sub_sess;
                extractData.xAxis(count).tc = x_Axis.sub(1).sess(1).tc;
                % Plot ICA time course
                plot(extractData.xAxis(count).tc, extractData.timecourses(count).tc, timeCourseColor, 'parent', meanAxesH);
                %tag{count} = ['Mean over all subjects and sessions'];
                % plot if Calculate Standard Error Mean exists
                if ~strcmpi(calculateSEM, 'no')
                    hold on;
                    % Calculate Standard Error mean
                    SEM_meanSubjectTc_up = meanAll_SEM + meanAll_sub_sess;
                    plot(x_Axis.sub(1).sess(1).tc, SEM_meanSubjectTc_up, SEM_Color, 'parent', meanAxesH);
                    hold on;
                    SEM_meanSubjectTc_down = meanAll_sub_sess - meanAll_SEM ;
                    plot(x_Axis.sub(1).sess(1).tc, SEM_meanSubjectTc_down, SEM_Color, 'parent', meanAxesH);
                    hold off;
                    clear SEM_meanSubjectTc_up;
                    clear SEM_meanSubjectTc_down;
                    %tag{count} = ['Mean over all subjects and sessions + SEM'];
                end
                title(tag{count}, 'parent', meanAxesH);
                axis(meanAxesH, 'tight');
                set(meanAxesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
                set(meanAxesH, 'tag', tag{count});
                % Store all positions
                initialPosition(count, :) = [axesXOrigin, axesYOrigin, axesWidth, axesHeight]; %get(meanAxesH, 'position');
                axesHandles(count) = meanAxesH;
                lastOrigin = axesXOrigin;
                lastYOrigin = axesYOrigin;
                axesXOrigin = axesXOrigin + axesWidth + xOffset;
            end
            % end for checking sessions
        end
        % end for loop over sessions
    end
    % end for mean display
end
% end for checking data sets

try
    close(helpHandle);
catch
end

drawnow;

% Plot sliders if possible
% Add Horizontal Slider if the x position of the axis is near 1
if initialPosition(end, 1) > (1 - axesWidth - xOffset)
    % Horizontal Slider Origin & dimensions
    sliderXOrigin = 0; sliderYOrigin = 0;
    horzSliderWidth = 0.97; horzSliderHeight = 0.03;
    horzSliderPos = [sliderXOrigin sliderYOrigin horzSliderWidth horzSliderHeight];
    % Store this information in userdata for both the sliders
    sliderData.numSubjects = numSubjects;
    sliderData.numSessions = numSessions;
    sliderData.initialPositions = initialPosition;
    sliderData.axesHandles = axesHandles;
    sliderData.tag = tag;
    maxVal = abs(horzSliderWidth - xOffset - axesWidth - initialPosition(end, 1));
    minVal = 0;
    % define slider steps
    horzSliderH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'slider', 'position', ...
        horzSliderPos, 'tag', 'horizontalSlider', 'callback', {@horzSliderCallback, graphicsHandle}, 'userdata', ...
        sliderData, 'min', minVal, 'max', maxVal, 'value', minVal, 'sliderstep', sliderStep);

end

% Add Vertical Slider if necessary (axes position is less than 0.01 and
% number of subjects is 1
if initialPosition(end, 2) < 0.01
    % Vertical Slider Origin & dimensions
    if exist('horzSliderHeight', 'var')
        verSliderXOrigin = sliderXOrigin + horzSliderWidth; verSliderYOrigin = horzSliderHeight;
    else
        verSliderXOrigin = 0.97; verSliderYOrigin = 0;
    end
    verSliderWidth = 1 - verSliderXOrigin; verSliderHeight = 1 - verSliderYOrigin;
    vertSliderPos = [verSliderXOrigin verSliderYOrigin verSliderWidth verSliderHeight];
    % Store this information in userdata for both the sliders
    sliderData.numSubjects = numSubjects;
    sliderData.numSessions = numSessions;
    sliderData.initialPositions = initialPosition;
    sliderData.axesHandles = axesHandles;
    sliderData.tag = tag;
    % Plot slider
    maxVal = 0;
    minVal = -abs(initialPosition(end, 2) - yOffset);
    vertSliderH = icatb_uicontrol('parent', graphicsHandle, 'units', 'normalized', 'style', 'slider', 'position', ...
        vertSliderPos, 'tag', 'verticalSlider', 'max', maxVal, 'min', minVal, 'value', maxVal, 'userdata', ...
        sliderData, 'callback', {@verSliderCallback, graphicsHandle}, 'sliderstep', sliderStep);
end

extractData.tag = tag;
set(extract_tc_Menu, 'callback', {@extractTc_Callback, graphicsHandle}, 'userdata', extractData);

set(graphicsHandle, 'WindowButtonDownFcn', @buttonDownFcn);

% Function callbacks

% Vertical Slider Callback
function verSliderCallback(handleObj, evd, handles)

set(handles, 'pointer', 'watch');

% Get the user data and value
sliderData = get(handleObj, 'userdata');
scrollValue = get(handleObj, 'value');

% Initial position and axes handles
initialPosition = sliderData.initialPositions;
axesHandles = sliderData.axesHandles;

% Update axes positions
figPos = get(axesHandles, 'position');
% Convert cell to MAT
figPos = cell2mat(figPos);
figPos(:, 2) = initialPosition(:, 2) - scrollValue;
% Convert figpos to cell
figPos = mat2cell(figPos, ones(size(figPos, 1), 1), 4);

% Set position to all handles
set(axesHandles, {'position'}, figPos);

drawnow;

set(handles, 'pointer', 'arrow');


% Horizontal Slider Callback
function horzSliderCallback(handleObj, event_data, handles)

set(handles, 'pointer', 'watch');

% Get the user data and value
sliderData = get(handleObj, 'userdata');
scrollValue = get(handleObj, 'value');

% Initial position and axes handles
initialPosition = sliderData.initialPositions;
axesHandles = sliderData.axesHandles;

% Update axes positions
figPos = get(axesHandles, 'position');
% Convert cell to MAT
figPos = cell2mat(figPos);
figPos(:, 1) = initialPosition(:, 1) - scrollValue;
% Convert figpos to cell
figPos = mat2cell(figPos, ones(size(figPos, 1), 1), 4);

% Set position to all handles
set(axesHandles, {'position'}, figPos);

drawnow;

set(handles, 'pointer', 'arrow');

function extractTc_Callback(hObject, event_data, handles)
% function to extract ICA time courses

% user data
userData = get(hObject, 'userdata');

% time courses and their title
timecourses = userData.timecourses;
tag = userData.tag;
xAxis = userData.xAxis;

figureData = repmat(struct('tag', '', 'yAxis', [], 'xAxis', []), length(userData.tag), 1);
% extract time courses
for ii = 1:length(userData.tag)
    figureData(ii).tag = tag{ii};
    figureData(ii).yAxis = timecourses(ii).tc;
    figureData(ii).xAxis = xAxis(ii).tc;
end

% open input dialog box
prompt = {'Enter a valid file Name:'};
dlg_title = 'Save timecourse data as';
num_lines = 1;
def = {'ICA_Timecourses'};
% save the file with the file name specified
fileName = icatb_inputdlg2(prompt, dlg_title, num_lines, def);

if ~isempty(fileName)
    [pathstr, fName, extn] = fileparts(fileName{1});
    if isempty(pathstr)
        pathstr = pwd;
    end
    extn = '.mat';
    location_file = fullfile(pathstr, [fName, extn]);
    icatb_save(location_file, 'figureData');
    disp(['File saved as ', location_file]);
end


function buttonDownFcn(hObject, event_data, handles)
% Execute button down fcn

set(hObject, 'pointer', 'watch');

% Identify selectionType
selectionType = get(hObject, 'SelectionType');

cH = gco;
objType = get(gco, 'Type');

if strcmpi(selectionType, 'normal') || strcmpi(selectionType, 'open')
    if ~strcmpi(objType, 'axes')
        if strcmpi(get(get(cH, 'parent'), 'Type'), 'axes')
            cH = get(cH, 'parent');
        else
            set(hObject, 'pointer', 'arrow');
            return;
        end
    end
    graphicsHandle = icatb_getGraphics('Timecourse', 'normal', 'time_course');
    new_handle = copyobj(cH, graphicsHandle);
    set(new_handle, 'position', [0.15 0.15 0.75 0.75]);
    set(graphicsHandle, 'name', get(get(new_handle, 'Title'), 'string'));
end

set(hObject, 'pointer', 'arrow');