function [varargout] = icatb_convertTo_struct(varargin)

% convert the concatenated time courses to structures depending upon the
% data sets and diffTimePoints vector
%
% Input:
% 1. timecourse - vector of length N total time points
% 2. modelTimecourse - matrix of dimensions N by M where M is the
% regressors
% 3. numSubjects - Number of subjects
% 4. numSessions - Number of sessions
% 5. diffTimePoints - vector containing the lengths of different time
% courses concatenated
% Output:
% 1. timecourseStruct - structure of the form
% timecourseStruct.sub().sess().tc
% 2. modelTimecourseStruct - structure of the form
% modelTimecourseStruct.sub().sess().tc
modelTimecourseStruct = [];
diffTimePoints = [];
modelTimecourse = [];
for ii = 1:2:nargin
    if strcmp(lower(varargin{ii}), 'timecourse')
        timecourse = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'model')
        modelTimecourse = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'numsubjects')
        numSubjects = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'numsessions')
        numSessions = varargin{ii + 1};
    elseif strcmp(lower(varargin{ii}), 'difftimepoints')
        diffTimePoints = varargin{ii + 1};
    end
end

%%%%%%%%%% Check for the existence of the vars %%%%%%%%%%%%%%
nPoints = length(timecourse);
if size(timecourse, 1) ~= nPoints
    timecourse = timecourse';
end

% time course vector
if prod(size(timecourse)) ~= nPoints
    error('Time course must be a vector');
end

if ~isempty(modelTimecourse)
    % model timecourse
    if size(modelTimecourse, 1) ~= nPoints
        modelTimecourse = modelTimecourse';
    end
end

% apply defaults to number of subjects and sessions
if ~exist('numSubjects', 'var')
    numSubjects = 1;
end
if ~exist('numSessions', 'var')
    numSessions = 1;
end

if isempty(diffTimePoints)
    nTimePoints = ceil(nPoints / (numSubjects*numSessions));
    diffTimePoints = repmat(nTimePoints, 1, numSubjects*numSessions);
end

% loop over structure
count = 0;
startTp = 1; % starting point
for ii = 1:numSubjects
    for jj = 1:numSessions
        count = count + 1;
        endTp = sum(diffTimePoints(1:count)); % current vector end point
        timecourseStruct.sub(ii).sess(jj).tc = timecourse(startTp:endTp, 1);
        if ~isempty(modelTimecourse)
            modelTimecourseStruct.sub(ii).sess(jj).tc = modelTimecourse(startTp:endTp, :);
        end
        x_Axis.sub(ii).sess(jj).tc = (1:length(timecourseStruct.sub(ii).sess(jj).tc))';
        startTp = endTp + 1;
    end
end
clear timecourse;
varargout{1} = timecourseStruct;
varargout{2} = modelTimecourseStruct;
varargout{3} = x_Axis;