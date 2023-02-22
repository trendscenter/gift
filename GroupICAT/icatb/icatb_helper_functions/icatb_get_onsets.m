function [data_sessionNumber, onset_number] = icatb_get_onsets(selectedRegressor, spmVar, data_sessionNumber)
% gets the session number and onset based on the regressor name
% 
% Input: 
% 1. selectedRegressor - name of the regressor
% 2. spmVar - SPM structure containing temporal information
% 3. flag_sessNum_onset - get session number or onset or both and should be
% one of these ('flag_sessNum_onset', 'sess_number', 'onset')
% 4. data_sessionNumber - optional parameter (in case of session specific
% regressors)
% 
% Output:
% data_sessionNumber and onset_number

% By default:
if ~exist('flag_sessNum_onset', 'var')
    flag_sessNum_onset = 'session_number&onset';
end

obtainedString = strread(selectedRegressor, '%s', 'delimiter', '('); % find the strings
obtainedString = strread(obtainedString{2}, '%s', 'delimiter', ')'); % get the second string
if ~exist('data_sessionNumber', 'var')
    data_sessionNumber = str2num(obtainedString{1}); % session number 
end    
stringObtained = strread(obtainedString{2}, '%s', 'delimiter', '*'); % get the strings 
ref_name_compare = deblank(stringObtained{1}); % final string

% check the session number
if length(spmVar.Sess) < data_sessionNumber
    data_sessionNumber = length(spmVar.Sess);
    %disp('Changing session number ...');
    %disp('Session number specified exceeds the structure SPM.Sess');
end

% get the onsets here
numOnsets = length(spmVar.Sess(data_sessionNumber).U); % number of onsets
onset_number = 1; % initialise onset number
% loop over number of onsets
for ii = 1:numOnsets
    getNames = spmVar.Sess(data_sessionNumber).U(ii).name;
    allNames = str2mat(getNames); % all the names
    matchedIndex = strmatch(lower(ref_name_compare), lower(allNames), 'exact');
    if ~isempty(matchedIndex)
        onset_number = ii;
        break;
    end
end
