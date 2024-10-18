clear all
tStart = tic;

% add requirements to path
addpath( '../bin/2019_03_03_BCT' )       % Brain connectivity toolbox

% GICA example with fbirn dataset
clear params;
params.param_file = '/data/users2/salman/projects/fBIRN/current/data/ICAresults_C100_fbirn/fbirnp3_rest_ica_parameter_info.mat';
params.outpath = '../results/fbirn_debug/';
params.fit_method = 'mnr';
params.n_corr = 3;
params.skip_noise = 0;
params.skip_anatomical = 0;
params.skip_functional = 0;
params.noise_training_set = 'pre_fbirn_sub';
params.anatomical_atlas = 'aal';
params.threshold = 3;
params.functional_atlas = 'yeo_buckner';
% params.functional_atlas = 'gordon2016';
% params.functional_atlas = 'caren';
disp( 'Running the autolabeller on FBIRN dataset' )
label_auto_main( params );

% Spatial map example with neuromark template
clear params;
params.sm_path = '/data/mialab/competition2019/NetworkTemplate/NetworkTemplate_High_VarNor.nii';
params.mask_path = '/data/mialab/competition2019/NetworkTemplate/Mask.img';
params.outpath = '../results/neuromark_debug/';
params.fit_method = 'mnr';
params.n_corr = 3;
params.skip_noise = 0;
params.skip_anatomical = 0;
params.skip_functional = 0;
params.noise_training_set = 'pre_aggregate';
params.anatomical_atlas = 'aal';
params.threshold = 3;
params.functional_atlas = 'yeo_buckner';
% params.functional_atlas = 'gordon2016';
% params.functional_atlas = 'caren';
disp( 'Running the autolabeller on NeuroMark dataset' )
label_auto_main( params );

tEnd = toc(tStart)
