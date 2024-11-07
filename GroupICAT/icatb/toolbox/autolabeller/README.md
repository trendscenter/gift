[![DOI](https://zenodo.org/badge/253869566.svg)](https://zenodo.org/badge/latestdoi/253869566)

# Autolabeller

This tool can automatically classify noisy spatial maps of brain activity, and generate anatomical and functional labels of the spatial maps and a reordered functional network connectivity matrix.

## Using the autolabeller

Example code can be found in `src/example_label_ic.m`. To run autolabeller, you need to add the requirements to your path:

    % add requirements to path
    addpath('/path/to/GIFTtoolbox')
    addpath('../bin/autolabeller/') % add the autolabeller folder
    
The following code can be used to run using ica_parameter.mat as an input file, containing all the necessary parameters of your ICA analysis.

    % group ICA example 
    clear params;
    params.param_file = '/path/to/ica_parameter_info.mat'; %input file
    params.outpath = '/results/path';
    params.fit_method = 'mnr';
    params.n_corr = 3; 
    params.skip_noise = 0;
    params.skip_anatomical = 0;
    params.skip_functional = 0;
    params.noise_training_set = 'pre_fbirn_sub';
    params.anatomical_atlas = 'aal';
    params.threshold = 3;
    params.functional_atlas = 'yeo_buckner';
    disp( 'Running the autolabeller on the selected dataset' )
    label_auto_main( params );

    % Spatial map example with the Neuromark template
    clear params;
    params.sm_path = '/path/to/SpatialMaps.nii';
    params.mask_path = './Mask.img';
    params.outpath = '/path/to/results/';
    params.fit_method = 'mnr';
    params.n_corr = 3;
    params.skip_noise = 0;
    params.skip_anatomical = 0;
    params.skip_functional = 0;
    params.noise_training_set = 'pre_aggregate';
    params.anatomical_atlas = 'aal';
    params.threshold = 3;
    params.functional_atlas = 'yeo_buckner';
    disp( 'Running the autolabeller on the selected dataset' )
    label_auto_main( params );

## Parameters & outputs

### Inputs
* `params.param_file` Location of the GICA parameter file
* `params.sm_path` Location of NIFTI data containing spatial maps. Use this if you are not running GICA.
* `params.outpath` Output directory
* `params.n_corr` paramter defininf how many ROI top correlations to calculate for anatomical/functional labeling. Default is set to 3
* `params.threshold` Threshold value for the spatial maps. Default is set to 3
* `params.skip_noise` Set to 0 by default. If you do not want to run or already completed artifact detection step, set to 1. 
* `params.skip_anatomical` Set to 0 by default. If you do not want to run or already completed anatomical labeling step, set to 1. 
* `params.skip_functional` Set to 0 by default. If you do not want to run or already completed functional labeling step, set to 1. 
* `params.noise_training_set` Choose dataset to use to train the noisecloud model. Options are:
    - `pre_fbirn_sub`: when both spatial maps and timecourses are available, as in a GIFT output
    - `pre_aggregate`: when only spatial maps are available
* `params.anatomical_atlas` Choose which atlas to use for anatomical labelling. Options: `aal`
* `params.functional_atlas` Choose which atlas to use for functional labelling. Options: `yeo_buckner`, `gordon2016`, `caren`. Default is `yeo_buckner`.

### Outputs
The following files are written into params.outpath folder:
* `network_labels.csv` is a vector of 0s or 1s corresponding to the input spatial maps; 0=artifact, 1=network
* `anatomical_labels.csv` containing the following columns:
    * `volume` 1-N where N is the number of input spatial maps
    * `network` a vector of 0s and 1s corresponding to the input spatial maps; 0=artifact, 1=network
    * `region_1`,`spatial_corr_1` AAL anatomical region with the highest spatial correlation to the spatial maps, and the correlation value
    * `region_2`,`spatial_corr_2`,`region_3`,`spatial_corr_3` AAL anatomical region with the second and third highest spatial correlations to the spatial maps, and the corresponding correlation values
* `functional_labels_[atlas_name].csv` has the following columns:
    * `volume` 1-N where N is the number of input spatial maps
    * `network` a vector of 0s and 1s corresponding to the input spatial maps; 0=artifact, 1=network
    * `region_1`,`spatial_corr_1` Functional parcellation from [atlas] with highest spatial correlation to the spatial maps, and the correlation value. Current available atlases are Yeo/BucknerLab, Gordon (2016), and CAREN 
    * `region_2`,`spatial_corr_2`,`region_3`,`spatial_corr_3` Functional parcellations with the second and third highest correlations to the spatial maps, and the corresponding correlation values
* `sorted_IC_idx_[atlas].csv` sorted index of the input spatial maps corresponding to the brain networks (artifact-related component indexes are removed)
* `sorted_fnc_[atlas].csv` sorted functional network connectivity (FNC) matrix of the brain networks 
* `nc` folder contains the noisecloud toolbox output (run under GIFT toolbox). It has the following files:
    * `*.nii` template `nii` files warped into the same space as the input spatial maps.
    * `nc_class_labels.txt` a vector of 0s and 1s corresponding to the input spatial maps; 0=artifact, 1=network
    * `training/testing_features.csv` contains the training/testing input data features used by the noisecloud toolbox in classification.  

## Result

The following figures are generated using the `example_plot_fnc.m` script. 
You can update the ICA parameter file and autolabeller output folder locations in the example codes above to generate new figures.
The script uses ICA parameter file to load the FNC from the ICA post-process result.

<img src="results/fbirn_nc_train_sub_th04/FBIRN_fnc_unsorted_yeo_buckner.png" alt="unsorted" width="150"/> <img src="results/fbirn_nc_train_sub_th04/FBIRN_fnc_reordered_yeo_buckner.png" alt="reordered" width="150"/> <img src="results/fbirn_nc_train_sub_th04/FBIRN_fnc_icn_yeo_buckner.png" alt="icn" width="150"/>

## Customizing the output

The autolabeler outputs can be easily updated based on visual observation as follows:
- Change the network(1)/noise(0) labels corresponding to the IC you want to update in `network_labels.csv`.
- Set `params.skip_noise = 1`
- Run the autolabeller with the original parameters again.

This will generate the updated anatomical/functional label files and IC order for the FNC matrix.

## Citation

Salman, M. S., Wager, T., Damaraju, E., Abrol, A., Vergara, V., Fu, Z., & Calhoun, V. (2021). An Approach to Automatically Label & Order Brain Activity/Component Maps. Brain Connectivity, brain.2020.0950. https://doi.org/10.1089/brain.2020.0950
