function icatb_ttest2_maps(files1, files2, outputDir, group1Name, group2Name, explicit_mask)
% Function to do two sample t-test on the selected files using SPM functions.
%
% Inputs:
% 1. files1 - Input files for group1
% 2. files2 - Input files for group2
% 3. outputDir - Output directory
% 4. group1Name - Group 1 name
% 5. group2Name - Group 2 name
% 6. explicit_mask - Explicit mask

if ~exist('outputDir', 'var')
    outputDir = pwd;
end

if ~exist('group1Name', 'var')
    group1Name = 'group 1';
end

if ~exist('group2Name', 'var')
    group2Name = 'group 2';
end

if (~exist('explicit_mask', 'var'))
    explicit_mask = '';
end

cd(outputDir);

spmMATFile = fullfile(outputDir, 'SPM.mat');

try
    if exist(spmMATFile, 'file')
        delete(spmMATFile);
    end
catch
end

load('icatb_spm_config_info.mat', 't2');

% Output directory
t2.jobs{1}.stats{1}.factorial_design.dir = {outputDir};

% Group1 files
t2.jobs{1}.stats{1}.factorial_design.des.t2.scans1 = cellstr(files1);

% Group2 files
t2.jobs{1}.stats{1}.factorial_design.des.t2.scans2 = cellstr(files2);

% Explicit mask
t2.jobs{1}.stats{1}.factorial_design.masking.em = cellstr(explicit_mask);

if (strcmpi(spm('ver'), 'spm8'))
    addpath(genpath(fileparts(which('spm.m'))));
    spm('Defaults', 'FMRI');
    jobs{1}.spm.stats.factorial_design = t2.jobs{1}.stats{1}.factorial_design;
else
    spm_defaults;
    global defaults;
    defaults.modality = 'FMRI';
    jobs = t2.jobs;
end

% Run job
spm_jobman('run', jobs);

% Load SPM.mat
spmInfo = load(spmMATFile);

% Estimate parameters
spmInfo.SPM = spm_spm(spmInfo.SPM);

% T-map contrast
spmInfo.SPM.xCon = spm_FcUtil('Set', [group1Name, ' - ', group2Name], 'T', 'c', [1;-1], spmInfo.SPM.xX.xKXs);

% SPM Contrasts
spm_contrasts(spmInfo.SPM);
