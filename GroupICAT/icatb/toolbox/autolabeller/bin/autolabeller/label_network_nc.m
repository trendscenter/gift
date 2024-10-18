% input: noisecloud training spatial map and timecourses
% output: 1 (network) or 0 (artifact)
function network_pred = label_network_nc( outpath, sesInfo, sm_path, tc_path, noise_training_set, threshold )
    disp('predicting networks')

    training_opts.pretrain = 0;
    training_opts.sm = [];
    training_opts.tc = [];

    if isa( noise_training_set, 'struct' )
        training_opts.sm = noise_training_set.sm; % Spatial maps
        training_opts.tc = noise_training_set.tc; % Timecourses
        % noise is lablled as 1 in noisecloud
        ic_meta = double( ~noise_training_set.class_labels );
    else
        switch noise_training_set
        case 'fbirn_sub'
            training_opts.sm = which('nc_training_sample_sm.nii'); % Spatial maps
            training_opts.tc = which('nc_training_sample_tc.nii'); % Timecourses
            % noise is lablled as 1 in noisecloud
            t1 = csvread( which( 'nc_training_labels.csv' ) );
            ic_meta = double( ~t1 );
        case 'pre_fbirn_sub'
            training_opts.pretrain = 1;
            % load precomputed features from file
            t1 = readtable( which('pre_fbirn_sub_th04.csv') );
            training_opts.feature_labels = t1.Properties.VariableNames;
            training_opts.features_norm = table2array( t1 );
            % use random FBIRN volume for registration
            training_opts.sm = which('fbirn_subxxx_component.nii');
            % noise is lablled as 1 in noisecloud
            t1 = csvread( which('pre_fbirn_sub_th04_labels.csv') );
            ic_meta = double( ~t1 );
        case 'pre_aggregate'
            training_opts.pretrain = 1;
            % load precomputed features from file
            t1 = readtable( which('pre_aggregate.csv') );
            training_opts.feature_labels = t1.Properties.VariableNames;
            training_opts.features_norm = table2array( t1 );
            % use HCP first volume for registration
            training_opts.sm = which('hcp_vol1_registration.nii');
            % noise is lablled as 1 in noisecloud
            t1 = csvread( which('pre_aggregate_labels.csv') );
            ic_meta = double( ~t1 );
        case 'debug'
            % load dummy 11 ICs from an FBIRN subject
            training_opts.sm = which('sm11_debug.nii'); % Spatial maps
            % noise is lablled as 1 in noisecloud
            ic_meta = ones( 1, 11 );
            t1 = [1 2 5 6 7 9 10];
            ic_meta( t1 ) = 0;
        otherwise
            % use fBIRN for training
            training_opts.sm = which('fbirnp3_rest_mean_component_ica_s_all_.nii'); % Spatial maps
            training_opts.tc = which('fbirnp3_rest_mean_timecourses_ica_s_all_.nii'); % Timecourses
            ic_meta_path = which( 'fBIRN_rsn.csv' );
            t1 = readtable( ic_meta_path );
            t1 = t1.IC;
            % noise is lablled as 1 in noisecloud
            ic_meta = ones( 1, 100 );
            ic_meta( t1 ) = 0;
        end
    end

    % other training params
    training_opts.regress_cov = [];
    training_opts.TR = 2;
    training_opts.class_labels = ic_meta; % Labels 
    
    if ~isempty( sesInfo )
        testing_opts.sm = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix '_mean_component_ica_s_all_.nii']);
        if ~isempty( training_opts.tc ) | training_opts.pretrain
            testing_opts.tc = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix '_mean_timecourses_ica_s_all_.nii']);
        else
            testing_opts.tc = [];
        end
    else
        testing_opts.sm = sm_path;
        if ~isempty( training_opts.tc ) | training_opts.pretrain
            testing_opts.tc = tc_path;
        else
            training_opts.tc = [];
            testing_opts.tc = [];
        end
    end

    testing_opts.regress_cov = [];
    testing_opts.TR = 2;

    nc_coregister = 1;
    if isempty( training_opts.sm )
        nc_coregister = 0;
    end

    [network_pred, fit_mdl, result_nc_classifier] = noisecloud_run(training_opts, testing_opts, 'convert_to_z', 'yes', 'outDir', outpath, 'coregister', nc_coregister, ...
        'iterations', 1, 'cross_validation', 10, 'threshold', threshold);
    
    % flip back because noise is lablled as 1 in noisecloud
    network_pred = ~network_pred;

    disp('done predicting network')
end

