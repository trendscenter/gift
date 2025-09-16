function icatb_ttest_maps(files, outputDir, group1Name)
% Function to do one sample t-test on the selected files using SPM functions.
%
% Inputs:
% 1. files - Input files
% 2. outputDir - Output directory
% 3. group1Name - Group Name
%

if ~exist('outputDir', 'var')
    outputDir = pwd;
end

if ~exist('group1Name', 'var')
    group1Name = 'group';
end

cd(outputDir);

spmMATFile = fullfile(outputDir, 'SPM.mat');

try
    if exist(spmMATFile, 'file')
        delete(spmMATFile);
    end
catch
end

load('icatb_spm_config_info.mat', 't1');

% Output directory
t1.jobs{1}.stats{1}.factorial_design.dir = {outputDir};

% Input files
t1.jobs{1}.stats{1}.factorial_design.des.t1.scans = cellstr(files);

if (strcmpi(spm('ver'), 'spm8'))
    addpath(genpath(fileparts(which('spm.m'))));
    spm('Defaults', 'FMRI');
    jobs{1}.spm.stats.factorial_design = t1.jobs{1}.stats{1}.factorial_design;
else
    spm_defaults;
    global defaults;
    defaults.modality = 'FMRI';
    jobs = t1.jobs;
end

% Run job
spm_jobman('run', jobs);

% Load SPM.mat
spmInfo = load(spmMATFile);

% Estimate parameters
spmInfo.SPM = spm_spm(spmInfo.SPM);

% T-map contrast
spmInfo.SPM.xCon = spm_FcUtil('Set', group1Name, 'T', 'c', 1, spmInfo.SPM.xX.xKXs);

% SPM Contrasts
spm_contrasts(spmInfo.SPM);
