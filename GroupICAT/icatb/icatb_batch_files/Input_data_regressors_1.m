% Batch file for running sorting GUI or option for entering the regressors through a file
% there are three cases for entering design matrix

% case 1: Same regressors applied over all data sets
% case 2: Different regressors over sessions but same over subjects
% case 3: Different regressors over subjects and sessions

% NOTE: All the regressors should be equal if selected for sessions or
% subjects


% case 1: the regressor names

regress_all = {'Sn(1) right*bf(1)', 'Sn(1) left*bf(1)'};

% case 2 and case 3 are as follows:
% case 2: Include only the regressors for subject 1 sessions.
% case 3: Include regressors depending upon the number of subjects and sessions.
regress_sub1_sess1 = {'Sn(1) right*bf(1)', 'Sn(1) left*bf(1)'};

regress_sub2_sess1 = {'Sn(1) right*bf(1)', 'Sn(1) left*bf(1)'};

regress_sub3_sess1 = {'Sn(1) left*bf(1)', 'Sn(1) right*bf(1)'};