function icatb_paired_t_test_maps(files1, files2, explicitMask, outputDir, contrastName)
%% Function to do paired t-test on the selected files using SPM functions.
%
% Inputs:
% 1. files1 - Input files for session 1 in a cell array
% 2. files2 - Input files for session 2 in a cell array
% 3. explicitMask - Explicit mask (optional)
% 4. outputDir - Output directory (optional)
% 5. contrastName - Contrast Name (optional)
%

files1 = cellstr(files1);
files2 = cellstr(files2);

if (~exist('explicitMask', 'var'))
    explicitMask = '';
end

if (~exist('outputDir', 'var'))
    outputDir = pwd;
end

if (~exist('contrastName', 'var'))
    contrastName = 'Condition 1 - Condition 2';
end

cd(outputDir);

spmMATFile = fullfile(outputDir, 'SPM.mat');

try
    if exist(spmMATFile, 'file')
        delete(spmMATFile);
    end
catch
end

load('icatb_spm_config_info.mat', 'pt');

% Output directory
pt.jobs{1}.stats{1}.factorial_design.dir = {outputDir};
pt.jobs{1}.stats{1}.factorial_design.masking.em = cellstr(explicitMask);
for nF = 1:length(files1)
    pt.jobs{1}.stats{1}.factorial_design.des.pt(1).pair(nF).scans = cellstr(char(deblank(files1{nF}), deblank(files2{nF})));
end

if (strcmpi(spm('ver'), 'spm8'))
    addpath(genpath(fileparts(which('spm.m'))));
    spm('Defaults', 'FMRI');
    jobs{1}.spm.stats.factorial_design = pt.jobs{1}.stats{1}.factorial_design;
else
    spm_defaults;
    global defaults;
    defaults.modality = 'FMRI';
    jobs = pt.jobs;
end

% Run job
spm_jobman('run', jobs);

% Load SPM.mat
spmInfo = load(spmMATFile);

% Estimate parameters
spmInfo.SPM = spm_spm(spmInfo.SPM);

names = char(spmInfo.SPM.xX.name);

ind1 = strmatch('Condition_{1}', names, 'exact');

ind2 = strmatch('Condition_{2}', names, 'exact');

c = zeros(size(names, 1), 1);
c(ind1) = 1;
c(ind2) = -1;

% T-map contrast
spmInfo.SPM.xCon = spm_FcUtil('Set', contrastName, 'T', 'c', c, spmInfo.SPM.xX.xKXs);

% SPM Contrasts
spm_contrasts(spmInfo.SPM);
