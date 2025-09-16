function icatb_formDesignMat(varargin)
% load data file and form a SPM.mat file
% data file can be any ascii file or .MAT file with variable containing the
% timecourse (m by n) where each column represents a time course
%
% Inputs:
% All inputs must be in pairs
% 1. fileName - Input file name (full file path).
% 2. nscan - Scans per session
% 3. regress_session - Regressors per session
% 4. outName - Output file name.
%
% Output: forms spm design matrix
%
% Input file contains time course variable where m denotes the sum of scans over sessions and n
% denotes sum of regressors over sessions.
%
% Example 1:
% fileName = which('visuomotor_regressors.dat'); nscan = 220; regress_session = 3;
% outName = 'visuo_SPM.mat';
%
% Example 2:
% fileName = which('aod_regressors.dat'); nscan = [249 249]; regress_session = [7 7];
% outName = 'aod_SPM.mat';
%
% Function call: icatb_formDesignMat('input_file_name', fileName,
% 'scans_per_session', nscan, 'regressors_per_session', regress_session,
% 'output_file_name', outName);
%

if nargin == 0
    % specify data file to load
    fileName = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'file', 'title', ...
        'Please select the regressor data', 'filter', '*.dat;*.mat');

    [timeCourse] = get_data_file(fileName);

    % open input dialog box
    prompt = {'Enter a valid file name:', 'Enter session time points as a row vector like [200 200]', ...
        'Regressors per session like [5 5]'};
    dlg_title = 'Design matrix information';
    num_lines = 1;
    def = {'SPM', num2str(size(timeCourse, 1)), num2str(size(timeCourse, 2))};

    % save the file with the file name specified
    getAnswers = icatb_inputdlg2(prompt, dlg_title, num_lines, def);

    if isempty(getAnswers)
        error('cannot proceed as the input dialog box is closed.');
    end

    % ouput file name
    outName = getAnswers{1};

    % get scans information
    nscan = str2num(getAnswers{2});

    % regressors per session
    regress_session = str2num(getAnswers{3});

else

    % loop over number of args
    for ii = 1:2:nargin
        if strcmpi(varargin{ii}, 'input_file_name')
            fileName = varargin{ii + 1}; % input file name
        elseif strcmpi(varargin{ii}, 'scans_per_session')
            nscan = varargin{ii + 1}; % Scans per session
        elseif strcmpi(varargin{ii}, 'regressors_per_session')
            regress_session = varargin{ii + 1}; % Regressors per session
        elseif strcmpi(varargin{ii}, 'output_file_name')
            outName = varargin{ii + 1}; % output file name
        end
        % end for checking args
    end
    % end loop over args

    if ~exist('fileName', 'var')
        error('file name is missing');
    end

    if ~exist('nscan', 'var')
        error('Scans per session variable is missing');
    end

    if ~exist('regress_session', 'var')
        error('Regressors per session variable is missing');
    end

    if ~exist('outName', 'var')
        error('Output file name is missing');
    end

    [timeCourse] = get_data_file(fileName);

end
% end for checking number of args

numSessions = length(nscan);

if numSessions ~= length(regress_session)
    error('Regressors per session must be of the same length as scans per session');
end

% get the current directory
[currentDir, fName, extn] = fileparts(outName);

% current directory
if isempty(currentDir)
    currentDir = fileparts(fileName);
    if isempty(currentDir)
        currentDir = pwd;
    end
end

% check the row length
if size(timeCourse, 1) ~= sum(nscan)
    error(['Row length (', num2str(size(timeCourse, 1)), ') of time course must be equal to the sum of the scans (' , ...
        num2str(sum(nscan)) , ') over sessions.']);
end

% check the even number of regressors
if size(timeCourse, 2) ~= sum(regress_session)
    error(['Column length of time course (', num2str(size(timeCourse, 2)), ') must be equal to the sum of the regressors over sessions (', ...
        num2str(sum(regress_session)), ')']);
end

count = 0;
% loop over sessions
for nSess = 1:numSessions
    % name regressors
    for ii = 1:regress_session(nSess)
        count = count + 1;
        % store it in SPM.Sess
        regressorName = ['Regressor ', num2str(ii)];
        % store it in SPM.xX.X
        snNames{count} = ['Sn(', num2str(nSess), ') ', regressorName];
        % form the SPM.Sess structure
        spmVar.Sess(nSess).U(ii).name = regressorName;
        % form the SPM onset
        spmVar.Sess(nSess).U(ii).ons = [1:1:nscan(nSess)];
    end
    % end for naming regressors

end
% end loop over sesions

% set the time course field
spmVar.xX.X = timeCourse;

% tr and units
spmVar.xY.RT = 1; spmVar.xBF.UNITS = 'scans';

% get the number of scans
spmVar.nscan = nscan;

% set the spm names
spmVar.xX.name = snNames;

% spm mat file
spmMatFile = fullfile(currentDir, [fName, '.mat']);


icatb_save(spmMatFile, 'spmVar');

disp(['Design matrix is created as ', spmMatFile]);


function [timeCourse] = get_data_file(fileName)
% get data from file

% file name
if isempty(fileName)
    error('regressor data file is not selected');
end

% load data
[getVar] = load(fileName);

% if the file is of type .MAT
if isstruct(getVar)
    % get field names
    field_names = fieldnames(getVar);
    timeCourse = getfield(getVar, field_names{1});
    clear getVar;
elseif isnumeric(getVar)
    % if the file is ascii
    timeCourse = getVar;
    clear getVar;
else
    error('Unknown data file format');
end
% end for checking


