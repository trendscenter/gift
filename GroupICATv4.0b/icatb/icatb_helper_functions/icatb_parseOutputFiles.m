function [subjectICAFiles, meanICAFiles, tmapICAFiles, meanALL_ICAFile, groupICAFiles, stdICAFiles] = ...
    icatb_parseOutputFiles(varargin)
% [subjectICAFiles, meanICAFiles, tmapICAFiles] = icatb_parseOutputFiles(icaOutputFiles)
% Input:
% icaOutputFiles(index).ses(sessionNumber).name(componentNumber,:);
%
% Output:
% meanICAFiles(sessionNumber).name(componentNumber,:)
% tmapICAFiles(sessionNumber).name(componentNumber,:)
% subjectICAFiles(subjectNumber).ses(sessionNubmer).name(componentNumber,:)
%

%load defaults
icatb_defaults;
global GROUP_ICA_INDEX;
global SUBJECT_ICA_INDEX;
global MEAN_INDEX;
global TMAP_INDEX;
global STD_INDEX;
global MEAN_ALL_INDEX;

meanICAFiles.name = [];
tmapICAFiles.name = [];
groupICAFiles.name = [];
stdICAFiles.name = [];
meanALL_ICAFile.name = [];

for i = 1:2:nargin
    if strcmp(lower(varargin{i}), 'numofsub')
        numberOfSubjects = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'numofsess')
        numberOfSessions = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'icaoutputfiles')
        icaOutputFiles = varargin{i + 1};
    elseif strcmp(lower(varargin{i}), 'flagtimepoints')
        flagTimePoints = varargin{i + 1};
    end
end

if ~exist('numberOfSubjects', 'var')
    error('Number of subjects doesn''t exist');
end

if ~exist('numberOfSessions', 'var')
    error('Number of sessions doesn''t exist');
end

if ~exist('flagTimePoints', 'var')
    flagTimePoints = 'same_time_points';
end


if(numberOfSubjects == 1 &  numberOfSessions == 1)
    
    % subject files
    subjectICAFiles(1).ses(1).name = icaOutputFiles(1).ses(1).name;
    
else
    if(numberOfSubjects > 1)
        %mean files
        for i=1:numberOfSessions
            meanICAFiles(i).name = icaOutputFiles(MEAN_INDEX).ses(i).name;
        end
        
        meanALL_ICAFile.name = icaOutputFiles(MEAN_ALL_INDEX).ses(1).name;
        
        %tmap files
        for i=1:numberOfSessions
            tmapICAFiles(i).name = icaOutputFiles(TMAP_INDEX).ses(i).name;
        end
        
        %std files
        for i=1:numberOfSessions
            stdICAFiles(i).name = icaOutputFiles(STD_INDEX).ses(i).name;
        end     
        
        %subject files
        for i=1:numberOfSessions
            for k=1:numberOfSubjects
                subjectICAFiles(k).ses(i).name = icaOutputFiles(k+(SUBJECT_ICA_INDEX-1)).ses(i).name;
            end
        end
        
    elseif(numberOfSubjects == 1 &  numberOfSessions > 1)
        meanALL_ICAFile.name = icaOutputFiles(1).ses(1).name;
        %subject files
        for i=1:numberOfSessions
            subjectICAFiles(1).ses(i).name = icaOutputFiles(2).ses(i).name;
        end
    end
    %end
    
end