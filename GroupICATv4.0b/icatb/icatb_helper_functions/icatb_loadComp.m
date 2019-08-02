function varargout = icatb_loadComp(sesInfo, compNumber, varargin)
%% Load selected subject components
%
% Inputs:
% 1. sesInfo - Parameter file or sesInfo data structure
% 2. compNumber - Component number or numbers.
% 3. varargin - Arguments passed in pairs name by value
%   a. detrend_no - Detrend level. Options are 0, 1, 2, and 3
%   b. subject_ica_files - Subject ICA Component files (image file names in
%   a data structure). Image files are used only the scaled MAT files don't
%   exist.
%   c. vars_to_load - Variables to load like {'ic', 'tc'}
%   d. average_runs - Average runs. Options are 0 and 1.
%   e. subjects - Subjects of interest.
%   f. sessions - Sessions of interest
%
% Outputs:
% Outputs depend on the variables of interest
%

if (ischar(sesInfo))
    pp = fileparts(sesInfo);
    if (isempty(pp))
        pp = pwd;
    end
    load(sesInfo);
    sesInfo.outputDir = pp;
end

zipContents.zipFiles = {};
zipContents.files_in_zip(1).name = {};
if ~isfield(sesInfo, 'zipContents')
    sesInfo.zipContents = zipContents;
end

conserve_disk_space = 0;
if (isfield(sesInfo, 'conserve_disk_space'))
    conserve_disk_space = sesInfo.conserve_disk_space;
end

averageRuns = 0;
subjects = (1:sesInfo.numOfSub);
sessions = (1:sesInfo.numOfSess);
detrendNumber = [];
truncateTp = 0;
covariates = '';
scansToInclude = [];

for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'detrend_no'))
        detrendNumber = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'subject_ica_files'))
        subjectICAFiles = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'vars_to_load'))
        vars_to_load = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'average_runs'))
        averageRuns = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'subjects'))
        subjects = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'sessions'))
        sessions = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'truncate_tp'))
        truncateTp = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'covariates'))
        covariates = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'scansToInclude'))
        scansToInclude = varargin{i + 1};
    end
end

% if (length(sessions) == 1)
%     averageRuns = 0;
% end


if (~isempty(covariates))
    covariates = cellstr(covariates);
end

if (~exist('compNumber', 'var'))
    compNumber = (1:sesInfo.numComp);
end

subjects = subjects(:)';
sessions = sessions(:)';
compNumber = compNumber(:)';

if (~exist('subjectICAFiles', 'var'))
    subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess);
end

modalityType = 'fmri';
if (isfield(sesInfo, 'modality'))
    modalityType = sesInfo.modality;
end

if (~strcmpi(modalityType, 'eeg'))
    compSetFields = {'ic', 'tc'};
else
    compSetFields = {'timecourse', 'topography'};
end

if (~exist('vars_to_load', 'var'))
    vars_to_load = compSetFields;
end

if (~iscell(vars_to_load))
    vars_to_load = cellstr(vars_to_load);
end

vars_to_load = intersect(vars_to_load, compSetFields);

if (isempty(vars_to_load))
    error('Wrong values passed to vars_to_load. Please check variable vars_to_load');
end

loadIC = ~isempty(strmatch(compSetFields{1}, vars_to_load, 'exact'));
loadTC = ~isempty(strmatch(compSetFields{2}, vars_to_load, 'exact'));

outDir = sesInfo.outputDir;

indices = zeros(1, length(subjects)*length(sessions));
endT = 0;
for nS = 1:length(subjects)
    startT = endT + 1;
    endT = endT + length(sessions);
    indices(startT:endT) = (subjects(nS) - 1).*sesInfo.numOfSess + sessions;
end

useCell = 0;
if (length(unique(sesInfo.diffTimePoints(indices))) ~= 1)
    useCell = 1;
end

useCell = (useCell && ~averageRuns && ~truncateTp);

numRuns = 1;
if (~averageRuns)
    numRuns = length(sessions);
end

tp = min(sesInfo.diffTimePoints(indices));

if (loadIC)
    if (useCell)
        IC = cell(length(subjects), length(sessions));
    else
        IC = zeros(length(subjects), numRuns, length(sesInfo.mask_ind), length(compNumber));
    end
end


if (loadTC)
    if (useCell)
        TC = cell(length(subjects), length(sessions));
    else
        TC = zeros(length(subjects), numRuns, tp, length(compNumber));
    end
end


countSub = 0;
for i = subjects
    
    countSub = countSub + 1;
    
    if (averageRuns)
        if (loadIC)
            tmpIC = zeros(length(sesInfo.mask_ind), length(compNumber));
        end
        if (loadTC)
            tmpTC = zeros(tp, length(compNumber));
        end
    end
    
    countSess = 0;
    for j = sessions
        countSess = countSess + 1;
        %% Load Components
        fileIn = fullfile(outDir, [sesInfo.calibrate_components_mat_file, num2str(i), '-', num2str(j), '.mat']);
        if (conserve_disk_space ~= 1)
            if (exist(fileIn, 'file'))
                [ic, tc] = loadMAT(fileIn, vars_to_load, compSetFields, compNumber, detrendNumber);
            else
                [ic, tc] = loadIm(sesInfo, subjectICAFiles(i).ses(j).name, loadTC, loadIC, compNumber, detrendNumber);
            end
        else
            [ic, tc] = loadIm(sesInfo, subjectICAFiles(i).ses(j).name, loadTC, loadIC, compNumber, detrendNumber);
        end
        
        if (~isempty(covariates))
            tc = regress_cov(tc, covariates{j + (i - 1)*sesInfo.numOfSess}, scansToInclude);
        end
        
        if (truncateTp)
            tc =  tc(1:tp, :);
        end
        
        if (useCell)
            %% Different timepoints
            if (loadIC)
                IC{countSub, countSess} = ic;
            end
            if (loadTC)
                TC{countSub, countSess} = tc;
            end
        else
            %% Same time points
            if (~averageRuns)
                if (loadIC)
                    IC(countSub, countSess, :, :) = ic;
                end
                if (loadTC)
                    TC(countSub, countSess, :, :) = tc;
                end
            else
                if (loadIC)
                    tmpIC = tmpIC + ic;
                end
                if (loadTC)
                    tmpTC = tmpTC + tc(1:tp, :);
                end
            end
        end
        clear ic tc;
    end
    
    %% Average runs
    if (averageRuns)
        if (loadIC)
            IC(countSub, 1, :, :) = (tmpIC/length(sessions));
        end
        if (loadTC)
            TC(countSub, 1, :, :) = (tmpTC/length(sessions));
        end
    end
    
    clear tmpIC tmpTC;
    
end

%% Remove singleton dimensions
if (loadIC && isnumeric(IC))
    IC = squeeze(IC);
end

if (loadTC && isnumeric(TC))
    TC = squeeze(TC);
end

if (loadIC && loadTC)
    varargout = {TC, IC};
else
    if (loadIC)
        varargout = {IC};
    else
        varargout = {TC};
    end
end


function [ic, tc] = loadIm(sesInfo, files, loadTC, loadIC, compNumber, detrendNumber)

currentFile = deblank(files(1, :));

zipFileName = {};
files_in_zip = {};
if ~exist(fullfile(sesInfo.outputDir, currentFile), 'file')
    [zipFileName, files_in_zip] = icatb_getViewingSet_zip(currentFile, [], 'real', sesInfo.zipContents);
    if (~isempty(zipFileName))
        icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), fullfile(sesInfo.outputDir, fileparts(currentFile)));
    end
end

files = icatb_fullFile('files', files, 'directory', sesInfo.outputDir);
files = icatb_rename_4d_file(files);
files = files(compNumber, :);

ic = [];
tc = [];

if (loadTC)
    tc = icatb_loadICATimeCourse(files);
    tc = tc(:, compNumber);
    if (~isempty(detrendNumber))
        tc = icatb_detrend(tc, 1, [], detrendNumber);
    end
end

if (loadIC)
    ic = icatb_read_data(files, [], sesInfo.mask_ind);
end

if (~isempty(zipFileName))
    icatb_delete_file_pattern(files_in_zip, sesInfo.outputDir);
end

function [ic, tc] = loadMAT(fileIn, vars_to_load, compSetFields, compNumber, detrendNumber)

ic = [];
tc = [];

info = load(fileIn, vars_to_load{:});
if (isfield(info, compSetFields{1}))
    ic = getfield(info, compSetFields{1});
    ic = ic(compNumber, :)';
end

if (isfield(info, compSetFields{2}))
    tc = getfield(info, compSetFields{2});
    tc = tc(:, compNumber);
    if (~isempty(detrendNumber))
        tc = icatb_detrend(tc, 1, [], detrendNumber);
    end
end

function tc = regress_cov(tc, file_name, scansToInclude)
%% Regress covariates from timecourses
%

file_name = deblank(file_name);

X = icatb_load_ascii_or_mat(file_name);

if (~exist('scansToInclude', 'var') || isempty(scansToInclude))
    scansToInclude = (1:size(X, 1));
end

scansToInclude(scansToInclude > size(X, 1)) = [];

if (isempty(scansToInclude))
    error(['Please check file numbers specified for file ', file_name]);
end

X = icatb_zscore(X);

X = X(scansToInclude, :);

if (size(X, 1) ~= size(tc, 1))
    error(['Please check the timepoints in file ', file_name]);
end

% Include temporal derivatives as well
%X = [X, [zeros(1, size(X, 2)); diff(X)]];

betas = pinv(X)*tc;

% Remove variance associated with the covariates
tc = tc - X*betas;
