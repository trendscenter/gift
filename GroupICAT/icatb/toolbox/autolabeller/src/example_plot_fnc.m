% sort the functional label files in order of domain names before running this for reproduction

clear all
tStart = tic;

addpath( genpath( '/trdapps/linux-x86_64/matlab/toolboxes/GroupICATv4.0b/' ) )      % GIFT toolbox
addpath( '../bin/my_icatb_plot_FNC.m' )

outpath = '../results/fbirn_nc_train_sub_th04/';
param_file = '/data/users2/salman/projects/fBIRN/current/data/ICAresults_C100_fbirn/fbirnp3_rest_ica_parameter_info.mat';
dataset_ = 'FBIRN';
% functional_atlas = 'yeo_buckner';
functional_atlas = 'gordon2016';
% functional_atlas = 'caren';

% outpath = '../results/cobre_nc_train_sub_th04/';
% param_file = '/data/users2/salman/projects/COBRE/current/results/ica_results_old/cobre1_ica_parameter_info.mat';
% dataset_ = 'COBRE';
% % functional_atlas = 'yeo_buckner';
% % functional_atlas = 'gordon2016';
% functional_atlas = 'caren';

fontname = 'Jost';

% load outputs
func_labels = readtable( fullfile( outpath, ['functional_labels_' functional_atlas '.csv'] ) );
func_labels = func_labels( func_labels.network==1, : );

sesInfo = load(param_file);
sesInfo = sesInfo.sesInfo;
num_IC = sesInfo.numComp;

% plot unsorted FNC for comparison
post_process = load( fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix '_postprocess_results.mat']) );
% fnc_unsorted = squeeze( mean( post_process.fnc_corrs_all ) );
t1 = size( post_process.fnc_corrs_all );
t2 = numel( t1( 1:end-2 ) );
fnc_unsorted = squeeze( mean( post_process.fnc_corrs_all, 1:t2 ) );

sorted_idx = func_labels.volume;
% plot FNC including noise
fnc = fnc_unsorted(sorted_idx, sorted_idx);

noise_idx = setdiff(1:num_IC, sorted_idx)';
full_idx = [sorted_idx; noise_idx];
max_fnc = max( abs( fnc_unsorted(:) ) );
% load module labels
[mod_names, t2, aff] = unique( func_labels.region_1, 'stable' );
mod_ = accumarray(aff, 1);
% add noise domain
mod_names = {'ICN', 'noise'};
mod_ = [length(sorted_idx); length(noise_idx)];

figure
my_icatb_plot_FNC(fnc_unsorted(full_idx, full_idx), [-max_fnc max_fnc], cell(1, num_IC), full_idx, gcf, 'Correlation', [], mod_, mod_names, 1);
title({['(E) ',dataset_,' dataset reordered FNC matrix']}, 'fontname', 'Jost', 'fontsize', 14)
set(gcf, 'color', 'w')
export_fig(fullfile(outpath, [dataset_ '_fnc_reordered_' functional_atlas '.png']), '-r150', '-p0.01')

% plot unsorted matrix
sorted_idx = 1:num_IC;
mod_ = num_IC;
mod_names = {''};

figure
my_icatb_plot_FNC(fnc_unsorted, [-max_fnc max_fnc], cell(1, num_IC), 1:num_IC, gcf, 'Correlation', [], mod_, mod_names, 1);
title({['(D) ',dataset_,' dataset unsorted FNC matrix']}, 'fontname', 'Jost', 'fontsize', 14)
ylabel('Components', 'fontweight', 'normal', 'fontsize', 14, 'fontname', fontname)
set(gcf, 'color', 'w')
export_fig(fullfile(outpath, [dataset_ '_fnc_unsorted_' functional_atlas '.png']), '-r150', '-p0.01')

% plot only ICN
sorted_idx = readmatrix( fullfile( outpath, ['sorted_network_idx_' functional_atlas '.csv'] ) );
max_fnc = max( abs( fnc_unsorted(:) ) );
% load module labels
[mod_names, t2, aff] = unique( func_labels.region_1, 'stable' );
mod_ = accumarray(aff, 1);

figure
my_icatb_plot_FNC(fnc_unsorted(sorted_idx, sorted_idx), [-max_fnc max_fnc], cell(1, num_IC), sorted_idx, gcf, 'Correlation', [], mod_, mod_names, 1);
title({['(A) ',dataset_,' dataset ICN FNC matrix'], ['Functional parcellation: ' strrep(functional_atlas, '_', '-')]}, 'fontname', 'Jost', 'fontsize', 14)
set(gcf, 'color', 'w')
export_fig(fullfile(outpath, [dataset_ '_fnc_icn_' functional_atlas '.png']), '-r150', '-p0.01')

close all

tEnd = toc(tStart)