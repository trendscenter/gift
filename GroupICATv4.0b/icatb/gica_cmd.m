function gica_cmd(varargin)
%% --GICA command line. Useful tool to run group ica with minimal options.
%
% Arguments must be passed in pairs
% Required:
%   1. --d or --data - Full file pattern or file name of the subject. You need to enter file names in single quotes. Otherwise matlab will ignore the spaces. For
%   multiple subjects use space to separate subjects like 'C:\sub01_vis\ns*img' 'C:\sub02_vis\ns*img'.
%
% Optional:
%   2. --o or --output - Output directory to place results. If omitted,
%   current directory is used.
%   3. --s or --sess - Number of sessions. Default is set to 1.
%   4. --preproc - Preproc type. Options are 1, 2, 3 or 4.
%       1. Remove Mean Per Timepoint
%       2. Remove Mean Per Voxel
%       3. Intensity Normalization
%       4. Variance Normalization
%   5. --m or --modality. Modality type. Options are fmri, eeg or smri.
%   6. --p or --prefix - Input prefix. Default is set to gica_cmd.
%   7. --n - Number of components extracted from the data. If
%   omitted, MDL tool will be used to estimate components from the data.
%   Mean of the estimated components is used as default.
%   8. --r or --recon - Backreconstruction type. Options are gica, str, gig-ica.
%   9. --a or --algorithm - ICA algorithm.
%   10. --pca - PCA algorithm. Options are:
%       1. Standard
%       2. Expectation Maximization
%       3. SVD
%       4. MPOWIT
%       5. STP
%   11. --icasso  - Stability analysis using ICASSO. Specify number of ica runs
%   followed by mode. Options for mode are randinit, bootstrap or both.
%   12. --mst - Stability analysis using MST. Specify number of ica runs.
%   13. --mask - Specify full file name of mask file or leave it empty for
%   default mask.
%   14. --dummy - Specify number of dummy scans or file numbers to include.
%   15. --parallel - Enter number of sessions/workers for doing analysis in
%   parallel.
%   16. --gtype - Group ica type. Options are spatial or temporal.
%   17. --reductions - Number of reduction steps. Max is 2.
%   18. --df - Enter number of components or degrees of freedom to retain in the first PCA step when 2 data reduction steps are used.
%   19. --performance - Options are 1, 2 and 3.
%       1 - Maximize Performance
%       2 - Less Memory Usage
%       3 - User Specified Settings
%   20. --templates - Enter fullfile path of template nifti file for constrained spatial ica
%   algorithms.
%   21. --display - Display results. Options are:
%       1 - HTML
%       2 - PDF
%
% COMMAND LINE USAGE:
%
% 1. File names are passed directly:
% gica_cmd --p visuo --o 'E:\test_GIFT\new_version\cmd' --data 'F:\Example Subjects\visuomotor\sub01_vis\ns*img' 'F:\Example Subjects\visuomotor\sub02_vis\ns*img' 'F:\Example Subjects\visuomotor\sub03_vis\ns*img' --n 30 --a infomax
%
% 2. Files are read through input text file:
% gica_cmd --p visuo --o 'E:\test_GIFT\new_version\cmd' --data files.txt --n 30 --a infomax


icatb_defaults;
global PREPROC_DEFAULT;
global BACKRECON_DEFAULT;
global PCA_DEFAULT;
global SCALE_DEFAULT;
global NUM_RUNS_GICA;

%% Parse inputs
inputs = varargin;
inds = icatb_good_cells(inputs);
inputs(inds==0)=[];
inds = icatb_good_cells(regexp(inputs, '--'));
matchInds = find(inds == 1);
vars = inputs(matchInds);


%% Initialise vars
if (isempty(PREPROC_DEFAULT))
    PREPROC_DEFAULT = 1;
end

if (isempty(PCA_DEFAULT))
    PCA_DEFAULT = 1;
end

if (isempty(BACKRECON_DEFAULT))
    BACKRECON_DEFAULT = 4;
end

if (isempty(SCALE_DEFAULT))
    SCALE_DEFAULT = 2;
end

outputDir = pwd;
numOfSess = 1;
preproc_type = PREPROC_DEFAULT;
modalityType = 'fmri';
prefix = 'gica_cmd';
backReconType = BACKRECON_DEFAULT;
algoType = 1;
pcaType = PCA_DEFAULT;
scaleType = SCALE_DEFAULT;
which_analysis = 1;
num_ica_runs = NUM_RUNS_GICA;
selMode = 'randinit';
maskFile = [];
dummy_scans = 0;
num_workers = 4;
analysisMode = 'serial';
group_ica_type = 'spatial';
numReductionSteps = 2;
perfType = 3; % user specified settings
display_results = 0;

for s = 1:length(matchInds)
    sInd = matchInds(s) + 1;
    sInd = min([sInd, length(inputs)]);
    try
        eInd = matchInds(s + 1) - 1;
    catch
        eInd =  length(inputs);
    end
    if (strcmpi(vars{s}, '--o') || strcmpi(vars{s}, '--output'))
        outputDir = inputs{sInd:eInd};
    elseif (strcmpi(vars{s}, '--d') || strcmpi(vars{s}, '--data'))
        filesP = inputs(sInd:eInd);
    elseif (strcmpi(vars{s}, '--s') || strcmpi(vars{s}, '--sess'))
        numOfSess = str2num(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--preproc'))
        tmp = inputs{sInd:eInd};
        tmp2 = [];
        try
            tmp2 = str2num(tmp);
        catch
        end
        if (~isempty(tmp2))
            preproc_type = tmp2;
        else
            preproc_type = tmp;
        end
    elseif (strcmpi(vars{s}, '--m') || strcmpi(vars{s}, '--modality'))
        modalityType = inputs{sInd:eInd};
    elseif (strcmpi(vars{s}, '--p') || strcmpi(vars{s}, '--prefix'))
        prefix = inputs{sInd:eInd};
    elseif (strcmpi(vars{s}, '--n'))
        numComp = str2num(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--r') || strcmpi(vars{s}, '--recon'))
        backReconType = inputs{sInd:eInd};
    elseif (strcmpi(vars{s}, '--a') || strcmpi(vars{s}, '--algorithm'))
        algoType = inputs{sInd:eInd};
    elseif (strcmpi(vars{s}, '--pca'))
        pcaType = inputs{sInd:eInd};
    elseif (strcmpi(vars{s}, '--icasso'))
        which_analysis = 2;
        tmp = inputs(sInd:eInd);
        num_ica_runs = str2num(tmp{1});
        try
            selMode = lower(tmp{2});
        catch
        end
        clear tmp;
    elseif (strcmpi(vars{s}, '--mst'))
        which_analysis = 3;
        num_ica_runs = str2num(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--mask'))
        maskFile = inputs{sInd:eInd};
    elseif (strcmpi(vars{s}, '--dummy'))
        dummy_scans = str2num(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--parallel'))
        analysisMode = 'parallel';
        num_workers = str2num(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--gtype'))
        group_ica_type = lower(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--df'))
        df = str2num(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--reductions'))
        numReductionSteps = str2num(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--performance'))
        tmp = inputs{sInd:eInd};
        tmp2 = [];
        try
            tmp2 = str2num(tmp);
        catch
        end
        if (~isempty(tmp2))
            perfType = tmp2;
        else
            perfType = tmp;
        end
    elseif (strcmpi(vars{s}, '--templates'))
        spatial_references = char(inputs{sInd:eInd});
    elseif (strcmpi(vars{s}, '--display'))
        display_results = str2num(inputs{sInd:eInd});
    end
end


disp('..............................................');
disp('............GICA COMMAND LINE ................');
disp('..............................................');
fprintf('\n');

if (ischar(algoType) && strcmpi(algoType, 'moo-icar'))
    algoType = 'gig-ica';
end

ica_algo = lower(cellstr(icatb_icaAlgorithm));
algoType = getIndex(algoType, ica_algo, 'ICA Algorithm');
%sesInfo.userInput.algorithm = algoType;

if (strcmpi(ica_algo{algoType}, 'semi-blind infomax'))
    error('GIFT cmd utility is not available for semi-blind ICA algorithms');
end

if (isempty(filesP))
    error('data parameter is missing. Command line usage is gica --o output_directory_name --data files1.nii --n 20 --a infomax');
end

if strcmpi(ica_algo{algoType}, 'gig-ica') || strcmpi(ica_algo{algoType}, 'constrained ica (spatial)')
    if (~exist('spatial_references', 'var'))
        error('Spatial references doesn''t exist for doing constrained ica');
    end
    inputData.refFiles = spatial_references;
end

[filesP, diffTimePoints] = listFiles(filesP);
numOfSub = (length(filesP)/numOfSess);
filesP = reshape(filesP, numOfSess, numOfSub)';

if (exist(outputDir, 'dir') ~= 7)
    mkdir(outputDir);
end

%% Fill inputData structure
inputData.inputFile = '';
inputData.modalityType = modalityType;
inputData.dataSelectionMethod = 4;
inputData.input_data_file_patterns = filesP;
inputData.dummy_scans = dummy_scans;
inputData.which_analysis = which_analysis;
inputData.outputDir = outputDir;
inputData.prefix = prefix;
inputData.maskFile = maskFile;
inputData.backReconType = backReconType;
inputData.preproc_type = preproc_type;

%% Group PCA performance settings
perfOptions = icatb_get_analysis_settings;
if (isnumeric(perfType))
    perfType = perfOptions{perfType};
end
perfType = lower(perfType);
inputData.perfType = perfType;

useTemporalICA = 0;
if (strcmpi(modalityType, 'fmri'))
    inputData.group_ica_type = lower(group_ica_type);
    useTemporalICA = strcmpi(group_ica_type, 'temporal');
end

if (numReductionSteps > 2)
    numReductionSteps = 2;
end

if (numOfSub*numOfSess == 1)
    numReductionSteps = 1;
end

if (useTemporalICA || strcmpi(ica_algo{algoType}, 'iva-gl') || strcmpi(ica_algo{algoType}, 'iva-l') || strcmpi(ica_algo{algoType}, 'gig-ica') || strcmpi(ica_algo{algoType}, 'constrained ica (spatial)'))
    numReductionSteps = 1;
end

inputData.numReductionSteps = numReductionSteps;

if (~exist('df', 'var'))
    df = min(diffTimePoints) - 1;
else
    df = min([df, min(diffTimePoints)]);
end

parallel_info.mode = lower(analysisMode);
parallel_info.num_workers = num_workers;
inputData.parallel_info = parallel_info;

doEstimation = 0;
if (~exist('numComp', 'var'))
    doEstimation = 1;
    if (numReductionSteps == 1)
        estimation_opts.PC1 = 'mean';
    else
        estimation_opts.PC1 = 'max';
        estimation_opts.PC2 = 'mean';
    end
    inputData.estimation_opts = estimation_opts;
else
    if (numReductionSteps == 1)
        inputData.numOfPC1 = numComp;
    else
        inputData.numOfPC1 = df;
        inputData.numOfPC2 = numComp;
    end
end

inputData.doEstimation = doEstimation;
inputData.pcaType = pcaType;
inputData.scaleType = scaleType;
inputData.algoType = algoType;

%% ICASSO Opts
if (which_analysis == 2)
    icasso_opts.num_ica_runs = max([2, num_ica_runs]);
    icasso_opts.sel_mode = selMode;
    icasso_opts.min_cluster_size = ceil(0.8*icasso_opts.num_ica_runs);
    icasso_opts.max_cluster_size = icasso_opts.num_ica_runs;
    inputData.icasso_opts = icasso_opts;
end

%% MST opts
if (which_analysis == 3)
    mst_opts.num_ica_runs = max([2, num_ica_runs]);
    inputData.mst_opts = mst_opts;
end

%% Setup analysis
param_file = icatb_read_batch_file(inputData);
load(param_file);

%% Run analysis
sesInfo = icatb_runAnalysis(sesInfo, 1);

%% Display results
if ((display_results ~= 0) && ~strcmpi(sesInfo.modality, 'eeg'))
    resultsF.formatName = display_results;
    icatb_report_generator(param_file, resultsF);
end


function [files, diffTimePoints] = listFiles(filesP)
%% list files

%% If text file is provided, read from file.
filesP = readTxtFile(filesP);

files = cell(length(filesP), 1);
diffTimePoints = zeros(1,  length(filesP));
for i = 1:length(filesP)
    fileP = filesP{i};
    [pathstr, fp, extn] = fileparts(fileP);
    fileContents = icatb_listFiles_inDir(pathstr, [fp, extn]);
    if (isempty(fileContents))
        error('Error:FilePattern', 'Please check file pattern %s as there are no files found\n', fileP);
    end
    fileListWithDir = icatb_fullFile('directory', pathstr, 'files', fileContents);
    fileListWithDir = icatb_rename_4d_file(fileListWithDir);
    files{i} = fileListWithDir;
    diffTimePoints(i) = size(fileListWithDir, 1);
end

function selOption = getIndex(selOption, options, titleStr)
%% Get selected index
%

if (ischar(selOption))
    ind = strmatch(lower(selOption), lower(cellstr(options)), 'exact');
    if (isempty(ind))
        error(['Unknown ', titleStr, '(', selOption, ') passed']);
    end
    selOption = ind(1);
end

function filesP = readTxtFile(filesP)

if (length(filesP) == 1)
    [dd, fN, extn] = fileparts(filesP{1});
    if (strcmpi(extn, '.txt'))
        fid = fopen(filesP{1}, 'r');
        if (fid == -1)
            error(['File ', filesP{1}, ' cannot be opened']);
        end
        filesP = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
        filesP = filesP{1};
        fclose(fid);
    end
end