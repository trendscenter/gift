% Batch file for running sorting GUI or option for entering the regressors through a file.

% There are three cases for entering regressors depending upon the way design matrix was selected
% in Display GUI.

% case 1: Same regressors applied over all data sets
% case 2: Different regressors over sessions but same over subjects
% case 3: Different regressors over subjects and sessions

% NOTE: All the regressors should be equal if selected for sessions or subjects.

% case 1: the regressor names
regress_all = {'Sn(1) TRG_PR_1_noslice*bf(1)', 'Sn(1) NOV_OM_1_noslice*bf(1)', 'Sn(1) STD_OM_1_noslice*bf(1)'};

% case 2 and case 3 are as follows:
% case 2: Include the regressors for subject 1 sessions.
% case 3: Include regressors depending upon the number of subjects and sessions.

% Subject 1 sessions
regress_sub1_sess1 = {'Sn(1) TRG_PR_1_noslice*bf(1)', 'Sn(1) NOV_OM_1_noslice*bf(1)'};
regress_sub1_sess2 = {'Sn(2) TRG_PR_2_noslice*bf(1)', 'Sn(2) NOV_OM_2_noslice*bf(1)'};

% Subject 2 sessions
regress_sub2_sess1 = {'Sn(1) TRG_PR_1_noslice*bf(1)', 'Sn(1) NOV_OM_1_noslice*bf(1)'};
regress_sub2_sess2 = {'Sn(2) TRG_PR_2_noslice*bf(1)', 'Sn(2) NOV_OM_2_noslice*bf(1)'};

% Subject 3 sessions
regress_sub3_sess1 = {'Sn(1) TRG_PR_1_noslice*bf(1)', 'Sn(1) NOV_OM_1_noslice*bf(1)'};
regress_sub3_sess2 = {'Sn(2) TRG_PR_2_noslice*bf(1)', 'Sn(2) NOV_OM_2_noslice*bf(1)'};