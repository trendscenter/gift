function icatb_open_subject_file
% Open subject file which contains information about data

% Form subject suffix
subSuffix = 'Subject';

% get subject file
sub_file = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
    ['*', subSuffix, '*.mat'], 'title', 'Select subject file');

% Load subject file
load(sub_file);


if ~exist('numOfSub', 'var') && ~exist('numOfSess', 'var') && ~exist('SPMFiles', 'var') && ~exist('files', 'var')
    % Check if the variables are present or not
    error(['Selected file: ', sub_file, ' is not a valid subject file']);
end


% Open setupICA analysis GUI
icatb_setup_analysis(sub_file);