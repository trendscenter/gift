function [out_file, comp_network_names] = icatb_merge_analyses(param_files, varargin)
%% Merge analyses and write only necessary files or symbolic links for doing mancova, dfnc or other group comparison tools.
% This tool works on nifti output only.
% Inputs:
%
% 1. param_files - Enter parameter files in a cell array
% 2. varargin - Arguments passed in pairs:
% a. merge_type:
% Option 1: 'stack_subjects' - Stack different spatially constrained ica analyses
% Option 2: 'stack_components' - Stack different model orders given the
% same subjects used across different analyses
% b. outputDir - Output Directory
% c. comp_network_names - Pass this variable if you are using
% 'stack_components' for 'merge_type'. Length of cell array must match
% number of ica parameter files.
% comp_network_names{1} = {'BG', [21, 30]; % Analysis 1 (30 components)
%                          'AUD', 17};
% comp_network_names{2} = {'SM', [10, 12]; % Analysis 2 (50 components)
%                          'DMN', [40, 42, 50]};
%
% Outputs:
% out_file - Output ica parameter file
% comp_network_names - Optional output. Component numbers are adjusted
% based on the final order of components in the out_file.
%

if (~iscell(param_files))
    param_files = cellstr(param_files);
end

% if (length(param_files) <= 1)
%     error('Need atleast two analyses');
% end


outDir = pwd;
merge_type = 'stack_subjects';

for n = 1:2:length(varargin)
    if (strcmpi(varargin{n}, 'merge_type'))
        merge_type = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'outputDir'))
        outDir = varargin{n + 1};
    elseif (strcmpi(varargin{n}, 'comp_network_names'))
        comp_network_names = varargin{n + 1};
    end
end

if (length(param_files) == 1)
    
    out_file = param_files{1};
    
    if (~exist('comp_network_names', 'var'))
        load(out_file);
        comp_network_names = {'ALL', (1:sesInfo.numComp)};
    end
    
    if (length(comp_network_names) == 1)
        comp_network_names = comp_network_names{1};
    end
    
    return;
    
end



numOfSess = 1;
numOfSub = 0;
numComp = 0;

ses = cell(length(param_files), 1);
tp = [];
files = [];
trs = [];

if (~exist('comp_network_names', 'var'))
    comp_network_names = cell(1, length(param_files));
end

all_comps = cell(1, length(param_files));

if (strcmpi(merge_type, 'stack_subjects'))
    missing_tr = zeros(1, length(param_files));
end

for n = 1:length(param_files)
    
    load(param_files{n});
    
    tmp_ses_outputDir = fileparts(param_files{n});
    if (isempty(tmp_ses_outputDir))
        tmp_ses_outputDir = pwd;
    end
    sesInfo.userInput.pwd = tmp_ses_outputDir;
    sesInfo.outputDir = tmp_ses_outputDir;
    
    if (n == 1)
        tmpComp = sesInfo.numComp;
        %tmpSub = sesInfo.numOfSub*sesInfo.numOfSess;
    end
    
    tmpSub = sesInfo.numOfSub*sesInfo.numOfSess;
    
    if (isempty(comp_network_names{n}))
        comp_network_names{n} = {'ALL', (1:sesInfo.numComp)};
    end
    
    if (strcmpi(merge_type, 'stack_subjects'))
        
        if isfield(sesInfo.userInput, 'TR')
            tmp_tr = sesInfo.userInput.TR;
            if (length(tmp_tr) == 1)
                tmp_tr = repmat(tmp_tr, 1, tmpSub);
            end
            trs = [trs, tmp_tr];
        else
            missing_tr(n) = 1;
        end
        
        numOfSub = tmpSub + numOfSub;
        numComp = sesInfo.numComp;
        if (tmpComp ~= (sesInfo.numComp))
            error('Number of components doesn''t match between analyses');
        end
        tp = [tp, sesInfo.diffTimePoints];
        start_tp = length(files) + 1;
        end_tp = start_tp + length(sesInfo.userInput.files) - 1;
        tmp_files = sesInfo.userInput.files;
        if (n == 1)
            files = tmp_files;
        else
            files(start_tp:end_tp) = tmp_files;
        end
        
    else
        
        files = sesInfo.userInput.files;
        if (tmpSub ~= (sesInfo.numOfSub*sesInfo.numOfSess))
            error('Number of data-sets doesn''t match between analyses');
        end
        numOfSub = sesInfo.numOfSub*sesInfo.numOfSess;
        
        tmp_cnames = comp_network_names{n};
        
        
        all_comps{n} = [tmp_cnames{:, 2}];
        
        %countCnums = numComp;
        for nca = 1:size(tmp_cnames, 1)
            tmp_cnames{nca, 1} = [tmp_cnames{nca, 1}, num2str(sesInfo.numComp)];
            tmp_cnums = tmp_cnames{nca, 2};
            comp_inds = numComp + (1:length(tmp_cnums));
            tmp_cnames{nca, 2} = comp_inds;
            numComp = numComp + length(tmp_cnums);
        end
        
        %numComp = numComp + length(all_comps{n});
        
        comp_network_names{n} = tmp_cnames;
        
        tp = sesInfo.diffTimePoints;
        
        
    end
    ses{n} = sesInfo;
    clear sesInfo;
    
end

sesInfo = ses{1};
sesInfo.userInput.files = files;
sesInfo.inputFiles = files;
sesInfo.userInput.numOfSub = numOfSub;
sesInfo.userInput.numOfSess = numOfSess;
sesInfo.userInput.numComp = numComp;
sesInfo.numComp = numComp;
sesInfo.numOfDataSets = sesInfo.userInput.numOfSub*sesInfo.userInput.numOfSess;
sesInfo.numOfSub = sesInfo.userInput.numOfSub;
sesInfo.numOfSess = sesInfo.userInput.numOfSess;
sesInfo.userInput.diffTimePoints = tp;
sesInfo.diffTimePoints = tp;
prefix = sesInfo.userInput.prefix;
prefix = [prefix, '_merge'];
sesInfo.userInput.prefix = prefix;
outputFiles = icatb_getOutputFileNames(sesInfo.userInput.prefix, numComp, numOfSub, numOfSess);
sesInfo.icaOutputFiles = outputFiles;

%[~, p, extn] = fileparts(sesInfo.userInput.param_file);
sesInfo.userInput.param_file = fullfile(outDir, [prefix, '_ica_parameter_info.mat']);
sesInfo.userInput.pwd = outDir;
sesInfo.outputDir = outDir;
sesInfo.param_files = param_files;


if (~isempty(trs))
    
    storeTr = 1;
    if (length(trs) ~= length(sesInfo.diffTimePoints))
        if (exist('missing_tr', 'var'))
            missingTRParamFiles = param_files(missing_tr == 1);
            for nMissing = 1:length(missingTRParamFiles)
                if (nMissing == 1)
                    disp('The following parameter files have missing TR ....');
                end
                disp(missingTRParamFiles{nMissing});
            end
            if (~isempty(missingTRParamFiles))
                storeTr = 0;
                disp('!!!!Concatenated TRs don''t match the number of data-sets. Not using TR information from the parameter files.');
            end
        end
        
    end
    
    if (storeTr)
        sesInfo.userInput.TR = trs;
        sesInfo.TR = trs;
    end
    
end

subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', numOfSub, 'numOfSess', numOfSess);

eT = 0;
allSubjectFiles = cell(numOfSub, length(ses));

for n = 1:length(ses)
    
    s1_param = ses{n};
    
    sfiles = icatb_parseOutputFiles('icaOutputFiles', s1_param.icaOutputFiles, 'numOfSub', s1_param.numOfSub, 'numOfSess', s1_param.numOfSess);
    
    tmpMask = zeros(s1_param.HInfo.DIM(1:3));
    tmpMask(s1_param.mask_ind) = 1;
    
    if (n == 1)
        mask = tmpMask;
    else
        mask = mask | tmpMask;
    end
    
    if (strcmpi(merge_type, 'stack_subjects'))
        %% Stack subjects
        
        for na = 1:s1_param.numOfSub
            for nb = 1:s1_param.numOfSess
                
                eT = eT + 1;
                tmpa = fullfile(outDir, subjectICAFiles(eT).ses(1).name);
                tmpb = fullfile(s1_param.outputDir, sfiles(na).ses(nb).name);
                
                tmpb = accumHdrFiles(tmpb);
                tmpa = accumHdrFiles(tmpa);
                
                for nLinks = 1:length(tmpb)
                    try
                        createlinks(tmpb{nLinks}, tmpa{nLinks});
                    catch
                    end
                end
                
            end
        end
        
        
    else
        %% STack components
        countKs = 0;
        for ksa = 1:length(sfiles)
            for ksb = 1:length(sfiles(ksa).ses)
                countKs = countKs + 1;
                tmp_files = icatb_rename_4d_file(icatb_fullFile('files', sfiles(ksa).ses(ksb).name, 'directory', ses{n}.outputDir));
                tmp_files = cellstr(tmp_files);
                allSubjectFiles{countKs, n} = tmp_files(all_comps{n});
            end
        end
        
        
    end
    
end


if (~strcmpi(merge_type, 'stack_subjects'))
    comp_network_names = cat(1, comp_network_names{:});
end


V = sesInfo.HInfo.V(1);
V.n(1) = 1;
V.fname = fullfile(outDir, [prefix, 'Mask.img']);
icatb_write_vol(V, double(mask));
mask_ind = find(mask == 1);
sesInfo.userInput.mask_ind = mask_ind;
sesInfo.mask_ind = mask_ind;
sesInfo.comp_network_names = comp_network_names;

%% Stack components and write out new component files
if (strcmpi(merge_type, 'stack_components'))
    
    for n = 1:length(subjectICAFiles)
        
        tmp_comp_names = subjectICAFiles(n).ses(1).name;
        tc = cell(1, size(allSubjectFiles, 2));
        ic = cell(1, size(allSubjectFiles, 2));
        for ncols = 1:size(allSubjectFiles, 2)
            tmpDat = icatb_read_data(allSubjectFiles{n, ncols}, [], mask_ind);
            tmpTC = icatb_loadICATimeCourse(char(allSubjectFiles{n, ncols}));
            tmpTC = tmpTC(:, all_comps{ncols});
            ic{ncols} = tmpDat';
            tc{ncols} = tmpTC;
        end
        
        ic = cat(1, ic{:});
        tc = cat(2, tc{:});
        
        icatb_saveICAData(tmp_comp_names, ic, tc, mask_ind, size(tc, 2), sesInfo.HInfo, 'real', [], outDir);
        
    end
    
end


fname = fullfile(outDir, [prefix, 'Subject.mat']);
files = sesInfo.inputFiles;
numOfSub = sesInfo.numOfSub;
numOfSess = sesInfo.numOfSess;

save(fname, 'numOfSub', 'numOfSess', 'files');
save(sesInfo.userInput.param_file, 'sesInfo');

out_file = sesInfo.userInput.param_file;

chkComponents = dir(fullfile(outDir, '*sub*comp*nii'));
chkTimecourses = dir(fullfile(outDir, '*sub*time*nii'));

if (isempty(chkTimecourses))
    chkTimecourses = dir(fullfile(outDir, '*sub*time*img'));
end

if (isempty(chkTimecourses))
    error('Subject component timecourses doesn''t exist. Cannot proceed with the analysis.');
end

if (isempty(chkComponents))
    chkComponents = dir(fullfile(outDir, '*sub*comp*img'));
end

try
    
	% Run group stats
    icatb_runAnalysis(sesInfo, 7);
	
catch
	
end

% if (~isempty(chkComponents))
%     % Run group stats
%     icatb_runAnalysis(sesInfo, 7);
% end


function createlinks(tmpb, tmpa)
%% Create symmbolic links
%

if ispc
    commandStr = ['mklink "', tmpa, '" ', tmpb];
else
    commandStr = ['ln -s "', tmpb, '" ', tmpa];
end

[status, message] = system(commandStr);

if (status)
    %error(message);
    disp(['Writing file ', deblank(tmpa(1, :)),  '...']);
    [dat, HInfo]=icatb_loadData(tmpb);
    icatb_write_nifti_data(tmpa, HInfo.V, dat);
end

%!mklink test_link.m Input_spatial_ica_test.m
%ln -s source_file symbolic_link


function files = accumHdrFiles(tmpb)
%% Get timecourse files and if any header files path

ta = getTimeCourseNaming(tmpb);

[~, fn, extn] = fileparts(deblank(tmpb(1, :)));
if (strcmpi(extn, '.img'))
    files = cell(size(tmpb, 1)*2 + 2, 1);
    count = 0;
    for nF = 1:size(tmpb, 1)
        
        currentFile = deblank(tmpb(nF, :));
        if (strcmp(extn, '.hdr'))
            tmpb_hdr = strrep(currentFile, extn, '.hdr');
            %ta_hdr = strrep(ta, extn, '.hdr');
        else
            tmpb_hdr = strrep(currentFile, extn, '.HDR');
            %ta_hdr = strrep(ta, extn, '.HDR');
        end
        
        count = count + 1;
        
        files{count} = currentFile;
        
        count = count + 1;
        files{count} = tmpb_hdr;
        
    end
    
    if (strcmp(extn, '.hdr'))
        ta_hdr = strrep(ta, extn, '.hdr');
    else
        ta_hdr = strrep(ta, extn, '.HDR');
    end
    
    count = count + 1;
    files{count} = ta;
    count = count + 1;
    files{count + 1} = ta_hdr;
    
    %files = {tmpb, tmpb_hdr, ta, ta_hdr};
else
    files = {tmpb, ta};
end



function timecourse_name = getTimeCourseNaming(P)

icatb_defaults;

global COMPONENT_NAMING;
global TIMECOURSE_NAMING;

P = icatb_parseExtn(deblank(P(1, :)));
% image extension
[pathstr_comp, file_comp, imExtn] = fileparts(P);
lastUnderScore = icatb_findstr(deblank(P(1,:)), '_');
lastUnderScore = lastUnderScore(end);
component_name = deblank(P(1, 1:lastUnderScore));
timecourse_name = strrep(component_name, COMPONENT_NAMING, TIMECOURSE_NAMING);
timecourse_name = [timecourse_name, imExtn];