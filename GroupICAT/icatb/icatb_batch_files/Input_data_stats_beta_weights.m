% Average beta weights across runs or sessions
% Options are 1 and 0
%
% Note: If you average runs, number of sessions will be set to 1 and it will not
% read from groupInfo variable.
averageRuns = 1;

% Design criteria
% Options are:
% 1 - One sample t-test (One condition)
% 2 - Two sample t-test (One condition)
% 3 - One Way Anova (Groups) (Multiple conditions)
% 4 - One Way Anova (Regressors) (Multiple conditions)
% 5 - Two Way Anova (Groups, Regressors) (Multiple conditions)
% 6 - Multiple Regression
desCriteria = 3;

% Group Info
% Enter group information in rows as group name, subject numbers, session numbers
%
% Note: You can also use matlab commands like load('sz.asc') to substitute for
% subject numbers and session numbers
groupInfo = {'SZ', (131:261), (1:4); ...
             'Healthy', (1:130), (1:4)};


% Selected groups
selGroups = [1 2];

% Selected conditions
selConditions = [1 3];

% Enter regressor files (Age, test scores) for using multiple regression in
% a cell array as shown below
multi_regress_files = {'C:\test_stats_beta\Age.txt';...
                       'C:\test_stats_beta\Scores.txt'};


%%%%%%%% OPTION FOR ONE WAY ANOVA (GROUPS) %%%%%%%%%%%%%%
%
% This option will let you do operation on regressors like Targets - Novels or
% Targets_run1 - Novels_run1 + Targets_run2 - Novels_run2 and the resultant regressor will be treated as  
% one condition.
%
% Length of equation for regressors must match the number of conditions.
% Leave empty if you don't want to use this option.
eq_regressors = [1 -1]; 
%
%%%%%%%% END FOR OPTION FOR ONE WAY ANOVA (GROUPS) %%%%%%%%%%%%%%


%%%%%% ANOVA CONTRASTS %%%%%%%%
% Contrast info
% Each row is a contrast name as shown below:
% contrastNames = {'SZ Vs Healthy'; 'Targets Vs Novels'};

contrastNames = {};

% Contrast matrix for anova where number of rows equal to number of contrasts
% Examples:
% 1. One way anova (Groups):
% contrastMatrix = [g1 g2 g3 etc];
% 2. One way anova (Regressors)
% contrastMatrix = [c1 c2 c3 c4 etc];
% 3. Two way anova (Groups, Regressors]
% contrastMatrix = [g1 g2 g3 c1 c2 c3 c4 etc];
% 

%                 g1  g2 c1 c2
%contrastMatrix = [1  -1  0  0; ...
%                  0   0  1 -1];

contrastMatrix = [];

%%%%%% END FOR ANOVA CONTRASTS %%%%%%%%