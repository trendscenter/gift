clear all
tStart = tic;

addpath( genpath( '../bin/MCIv4/' ) )
addpath( genpath( '../bin/CanlabCore/' ) )
addpath( '/data/users2/xxxx/projects/funfc/src/' )

% % fbirn
% outpath = '../results/fbirn_trained_by_neuromark/';
% param_file = '/data/users2/salman/projects/fBIRN/current/data/ICAresults_C100_fbirn/fbirnp3_rest_ica_parameter_info.mat';
% anatomical_atlas = 'aal';
% functional_atlas = 'yeo_buckner';

% cobre
outpath = '../results/cobre_trained_by_neuromark/';
param_file = '/data/users2/salman/projects/COBRE/current/results/ica_results_old/cobre1_ica_parameter_info.mat';
anatomical_atlas = 'aal';
functional_atlas = 'yeo_buckner';

structFile = '../bin/MCIv4/ch2better_whitebg_aligned2EPI_V4.nii';

mkdir( fullfile(outpath, 'sm_fig', 'nii') );

% load outputs
network_labels = readmatrix( fullfile( outpath, 'network_labels.csv' ) );
func_labels = table2cell( readtable( fullfile( outpath, ['functional_labels_' functional_atlas '.csv'] ) ) );
anat_labels = table2cell( readtable( fullfile( outpath, ['anatomical_labels_' anatomical_atlas '.csv'] ) ) );

sesInfo = load(param_file);
sesInfo = sesInfo.sesInfo;

% write output spatial maps with labels
agg_map_path = fullfile(sesInfo.outputDir, [sesInfo.aggregate_components_an3_file '.nii']);
mask_file = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix 'Mask.img']);
sm_dat = fmri_data( agg_map_path, mask_file, 'noverbose' );
n_vols = size( sm_dat.dat, 2 );

disp('plot the correlations')
corrs_ = readmatrix( fullfile(outpath, 'anatomical_correlations.csv') );
figure
fontname = 'Jost';
imagesc(corrs_);
axis image
xlabel('Atlas regions', 'fontname', fontname)
ylabel('Components', 'fontname', fontname)
title('Correlation between the spatial maps and AAL atlas regions', 'fontname', fontname)
c = colorbar( 'fontname', fontname );
ylabel(c, 'Correlation', 'fontname', fontname)
export_fig( fullfile(outpath, 'anatomical_correlations.png'), '-r150' )

for jj = 1:n_vols
    disp( ['plotting IC ' num2str(jj)] )
    % create title
    title_rsn = 'NO';
    title_anat = '';
    title_func = '';
    idx = find( [anat_labels{:,1}] == jj );
    if ( network_labels(jj) == 1 )
        title_rsn = 'YES'; 
    end
    % assume top 3 correlations 
    title_anat = [anat_labels{idx,3} ' (' num2str(anat_labels{idx,4}) ')']; 
    title_func = [func_labels{idx,3} ' (' num2str(func_labels{idx,4}) ')']; 
    
    title_ = ['RSN: ' title_rsn ';ANAT: ' title_anat ';FUNC: ' title_func];
    
    params = struct( ...
        'disable', 0, ...
        'data', sm_dat.dat(:, jj), ...
        'sesInfo', sesInfo, ...
        'structFile', structFile, ...
        'title', title_, ...
        'savefig', 1, ...
        'outpath', fullfile( outpath, 'sm_fig' ), ...
        'outname', ['fig' num2str(jj)] );
    funfc_mciv4_save(params);
    % break
end

close all

tEnd = toc(tStart)
