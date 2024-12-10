% Autolabeller function
%   Given a set of spatial maps (and optional timecourses from GICA output), 
%   autolabeller separates networks from artifacts (using the Noisecloud 
%   toolbox [1]), provides anatomical [2] and functional labels [3-4] and sorts 
%   the maps based on the functional labels (using BCT toolbox) [5].
% 
% Example usage: 
%   See src/example_label_id.m for examples.
% 
% Inputs
%   params.param_file
%       Location of GICA parameter file
%   params.sm_path
%       Location of NIFTI data containing spatial maps. Use this if you are not 
%       running GICA.
%   params.outpath
%       Output directory
%   params.n_corr
%       How many ROI top correlations to calculate for anatomical/functional
%       labeling. Default = 3
%   params.threshold
%       Threshold value for the spatial maps. Default = 3
%   params.skip_noise
%       If you do not want to run or already ran artifact detection step, set to
%       1. Otherwise set to 0 by default.
%   params.skip_anatomical
%       If you do not want to run or already ran anatomical labeling step, set to
%       1. Otherwise set to 0 by default.
%   params.skip_functional
%       If you do not want to run or already ran functional labeling step, set to
%       1. Otherwise set to 0 by default.
%   params.noise_training_set
%       Options: 'pre_fbirn_sub', 'pre_aggregate'
%       Which dataset to use to train the noisecloud model
%       pre_fbirn_sub: when both spatial maps and timecourses are available, as in a GIFT output
%       pre_aggregate: when only spatial maps are available
%   params.anatomical_atlas
%       Options: 'aal'
%       Which atlas to use for anatomical labeling
%   params.functional_atlas
%       Options: 'yeo_buckner', 'gordon2016', 'caren'
%       Which atlas to use for functional labeling. Default = 'yeo_buckner'
% 
% Outputs: the following files are written into params.outpath folder:
%   network_labels.csv
%       Network labels vector (0=artifact, 1=network) and probability that the component/spatial map is a network
%   anatomical_labels.csv
%       AAL anatomical region with highest correlations
%   functional_labels.csv
%       Buckner functional parcellations with highest correlations
%   sorted_IC_idx.csv
%       sorted IC index
%   sorted_fnc.csv
%       sorted FNC matrix
% 
% References:
%   [1] V. Sochat, K. Supekar, J. Bustillo, V. Calhoun, J. A. Turner, and D. L. Rubin, “A Robust Classifier to Distinguish Noise from fMRI Independent Components,” PLoS One, vol. 9, no. 4, Apr. 2014, doi: 10.1371/journal.pone.0095493.
%   [2] N. Tzourio-Mazoyer et al., “Automated Anatomical Labeling of Activations in SPM Using a Macroscopic Anatomical Parcellation of the MNI MRI Single-Subject Brain,” NeuroImage, vol. 15, no. 1, pp. 273–289, Jan. 2002, doi: 10.1006/nimg.2001.0978.
%   [3] B. T. T. Yeo et al., “The organization of the human cerebral cortex estimated by intrinsic functional connectivity,” J. Neurophysiol., vol. 106, no. 3, pp. 1125–1165, Sep. 2011, doi: 10.1152/jn.00338.2011.
%   [4] R. L. Buckner, F. M. Krienen, A. Castellanos, J. C. Diaz, and B. T. T. Yeo, “The organization of the human cerebellum estimated by intrinsic functional connectivity,” Journal of Neurophysiology, vol. 106, no. 5, pp. 2322–2345, Jul. 2011, doi: 10.1152/jn.00339.2011.
%   [5] M. Rubinov and O. Sporns, “Complex network measures of brain connectivity: Uses and interpretations,” NeuroImage, vol. 52, no. 3, pp. 1059–1069, Sep. 2010, doi: 10.1016/j.neuroimage.2009.10.003.
% 
% Citation:
%   https://doi.org/10.1101/2020.08.31.275578 


function label_auto_main( params )
    % add paths
    src_dir = fileparts( which('label_auto_main') );
    data_dir = fullfile( src_dir, 'data' );
    bin_dir = fullfile( src_dir, 'bin' );
    addpath( genpath( data_dir ) )
    addpath( genpath( bin_dir ) )
    
    % create output directory
    if ~exist( params.outpath, 'dir' )
        mkdir( fullfile( params.outpath ) )
        mkdir( fullfile( params.outpath, 'nc' ) )
    end
    if ~exist( [params.outpath filesep 'nc'], 'dir' )
        mkdir( fullfile( params.outpath, 'nc' ) )
    end

    % sanitize input
    if isfield( params, 'param_file' ) && ~isempty( params.param_file )
        % load ICA session info
        sesInfo = load( params.param_file );
        sesInfo = sesInfo.sesInfo;
        % load GICA post_process result 
        post_process = load( fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix '_postprocess_results.mat']) );
        % IC aggregate map path
        sm_path = fullfile(sesInfo.outputDir, [sesInfo.aggregate_components_an3_file '.nii']);
        tc_path = '';
        mask_path = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix 'Mask.img']);
        if ~exist( mask_path )
            mask_path = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix 'Mask.nii']);
        end
        flag_sort_fnc = 1;
        
        % for some reason autolabeller demands the serial postprocessing
        % directory so we create it 12/8/24
        s_dir_serial_postproc = fullfile(sesInfo.outputDir, [sesInfo.userInput.prefix '_postprocess_results']);
        if ~exist(s_dir_serial_postproc, 'dir')
            % create the contents in the dir as well 
            mkdir(s_dir_serial_postproc);
            icatb_postprocess_timecourses(sesInfo);
            % add a note
            s_dir_serial_postproc_readme = [s_dir_serial_postproc filesep 'readme_trends.txt'];
            [stat] = system(['echo ' '''serial postproc files created for autolabeller''' ' >  ' s_dir_serial_postproc_readme]);
        end
    else
        sesInfo = [];
        sm_path = params.sm_path;
        tc_path = '';
        if isfield( params, 'tc_path' )
            tc_path = params.tc_path;
        end
        mask_path = params.mask_path;
        flag_sort_fnc = 0;
    end

    if isfield( params, 'threshold' ) && ~params.threshold
        disp('Forcing threshold=1')
        params.threshold = 1;
    end

    if ~isfield( params, 'noise_training_set' )
        params.noise_training_set = 'fbirn';
    end

    % predict network labels (0=artifact, 1=network)
    if params.skip_noise
        nc_out = fullfile( params.outpath, 'network_labels.csv' );
        if exist( nc_out ) == 2
            network_labels = readmatrix( nc_out );
        else
            % For some reason user chose not to run noisecloud, so assume all components are resting-state networks.
            network_labels = 1;
        end
    else
        network_labels = label_network_nc( fullfile( params.outpath, 'nc' ), sesInfo, sm_path, tc_path, params.noise_training_set, params.threshold );
        csvwrite( fullfile(params.outpath, 'network_labels.csv'), network_labels )
    end
    
    % predict anatomical labels
    if params.skip_anatomical
        anat_labels = readtable( fullfile( params.outpath, ['anatomical_labels_' params.anatomical_atlas '.csv'] ) );
        anat_labels = [anat_labels.Properties.VariableNames; table2cell( anat_labels )];
    else
        [anat_labels, corrs_] = label_anatomical( sm_path, mask_path, params.threshold, network_labels, params.anatomical_atlas, params.n_corr, params.outpath );
        csvwrite( fullfile(params.outpath, 'anatomical_correlations.csv'), corrs_ )
    end
    
    % predict functional labels
    if params.skip_functional
        % check if the func label file exists
        fl = fullfile( params.outpath, ['functional_labels_' params.functional_atlas '.csv'] );
        if exist( fl )        
            func_labels = readtable( fl );
            func_labels = [func_labels.Properties.VariableNames; table2cell( func_labels )];
        end
    else
        func_labels = label_functional( sm_path, mask_path, params.threshold, network_labels, params.functional_atlas, params.n_corr, params.outpath );
    end

    % sort FNC
    if flag_sort_fnc
        % load unsorted FNC
        % works for session=1 only
        % fnc = squeeze( mean( post_process.fnc_corrs_all ) );
        if ~isempty(sesInfo)
	        s_dir_giftout = fileparts( params.param_file );
            clear post_process.fnc_corrs_all;
	        post_process.fnc_corrs_all=zeros(sesInfo.numOfSub,1,sesInfo.numComp,sesInfo.numComp);
	        for i_sub = 1:sesInfo.numOfSub	
	            tmp2=load([s_dir_giftout filesep sesInfo.userInput.prefix '_postprocess_results' filesep sesInfo.userInput.prefix '_post_process_sub_' sprintf( '%03d', i_sub )  '.mat'], 'fnc_corrs');
	            post_process.fnc_corrs_all(i_sub,1,:,:)=squeeze(tmp2.fnc_corrs);
                clear tmp2;
            end
	        % sort
            t1 = size( post_process.fnc_corrs_all );
            t2 = numel( t1( 1:end-2 ) );
            fnc = squeeze( mean( post_process.fnc_corrs_all, 1:t2 ) );
            [sorted_idx, network_fnc, order_] = sort_fnc( fnc, func_labels(2:end,1:3) );
    
            % sort the other labels
            if ~params.skip_anatomical
                anat_labels(2:end, :) = anat_labels( order_+1, : );
            end
            if ~params.skip_functional
                func_labels(2:end, :) = func_labels( order_+1, : );
            end
        else
	        disp('label_auto_main: Sorting needs the parameter file, which is missed');
        end
    end

    % write output
    if ~params.skip_anatomical
        writecell( anat_labels, fullfile(params.outpath, ['anatomical_labels_' params.anatomical_atlas '.csv']) )
    end
    if ~params.skip_functional
        writecell( func_labels, fullfile(params.outpath, ['functional_labels_' params.functional_atlas '.csv']) )
    end
    if flag_sort_fnc
        csvwrite( fullfile(params.outpath, ['sorted_network_idx_' params.functional_atlas '.csv']), sorted_idx )
        csvwrite( fullfile(params.outpath, ['sorted_fnc_' params.functional_atlas '.csv']), network_fnc )
    end

    

