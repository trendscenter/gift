[![DOI](https://zenodo.org/badge/253869566.svg)](https://zenodo.org/badge/latestdoi/253869566)

# Autolabeller

This tool can automatically classify noisy spatial maps of brain activity, and generate anatomical and functional labels of the spatial maps and a reordered functional network connectivity matrix.

## Prerequisites

Autolabeller is written in Matlab™ and requires several Matlab toolboxes to run. Please download the following toolboxes and add to your Matlab path.

- [BCT Toolbox](https://sites.google.com/site/bctnet/) (March 2019 release)

## Using the autolabeller

Example code can be found in `src/example_label_ic.m`.

    % add requirements to path
    addpath( genpath( '../bin/CanlabCore' ) )       % Canlab toolbox
    addpath( '../bin/2019_03_03_BCT' )       % Brain connectivity toolbox

    % GICA example with fbirn dataset
    clear params;
    params.param_file = './fbirnp3_rest_ica_parameter_info.mat';
    params.outpath = './results/fbirn/';
    params.fit_method = 'mnr';
    params.n_corr = 3;
    params.skip_noise = 0;
    params.skip_anatomical = 0;
    params.skip_functional = 0;
    params.noise_training_set = 'pre_fbirn_sub';
    params.anatomical_atlas = 'aal';
    params.threshold = 3;
    params.functional_atlas = 'yeo_buckner';
    disp( 'Running the autolabeller on FBIRN dataset' )
    label_auto_main( params );

    % Spatial map example with the Neuromark template
    clear params;
    params.sm_path = './NetworkTemplate_High_VarNor.nii';
    params.mask_path = './Mask.img';
    params.outpath = './results/neuromark/';
    params.fit_method = 'mnr';
    params.n_corr = 3;
    params.skip_noise = 0;
    params.skip_anatomical = 0;
    params.skip_functional = 0;
    params.noise_training_set = 'pre_aggregate';
    params.anatomical_atlas = 'aal';
    params.threshold = 3;
    params.functional_atlas = 'yeo_buckner';
    disp( 'Running the autolabeller on NeuroMark dataset' )
    label_auto_main( params );

## Parameters & outputs

### Inputs
* `params.param_file` Location of GICA parameter file
* `params.sm_path` Location of NIFTI data containing spatial maps. Use this if you are not running GICA.
* `params.outpath` Output directory
* `params.n_corr` How many ROI top correlations to calculate for anatomical/functional labeling. Default = 3
* `params.threshold` Threshold value for the spatial maps. Default = 3
* `params.skip_noise` If you do not want to run or already ran artifact detection step, set to 1. Otherwise set to 0 by default.
* `params.skip_anatomical` If you do not want to run or already ran anatomical labeling step, set to 1. Otherwise set to 0 by default.
* `params.skip_functional` If you do not want to run or already ran functional labeling step, set to 1. Otherwise set to 0 by default.
* `params.noise_training_set` Which dataset to use to train the noisecloud model. Options: `pre_fbirn_sub`, `pre_aggregate`
    - `pre_fbirn_sub`: when both spatial maps and timecourses are available, as in a GIFT output
    - `pre_aggregate`: when only spatial maps are available
* `params.anatomical_atlas` Which atlas to use for anatomical labeling. Options: `aal`
* `params.functional_atlas` Which atlas to use for functional labeling. Options: `yeo_buckner`, `gordon2016`, `caren`. Default = `yeo_buckner`.

### Outputs
The following files are written into params.outpath folder:
* `network_labels.csv` is a vector of 0/1 corresponding to the input spatial maps; 0=artifact, 1=network
* `anatomical_labels.csv` has the following columns:
    * `volume` 1-N where N is the number of input spatial maps
    * `network` a vector of 0/1 corresponding to the input spatial maps; 0=artifact, 1=network
    * `region_1`,`spatial_corr_1` AAL anatomical region with the highest spatial correlation to the spatial maps, and the correlation value
    * `region_2`,`spatial_corr_2`,`region_3`,`spatial_corr_3` AAL anatomical region with the second and third highest spatial correlations to the spatial maps, and the corresponding correlation values
* `functional_labels_[atlas].csv` has the following columns:
    * `volume` 1-N where N is the number of input spatial maps
    * `network` a vector of 0/1 corresponding to the input spatial maps; 0=artifact, 1=network
    * `region_1`,`spatial_corr_1` Functional parcellation from [atlas] with highest spatial correlation to the spatial maps, and the correlation value. Current available atlas are Yeo/BucknerLab, Gordon (2016), and CAREN 
    * `region_2`,`spatial_corr_2`,`region_3`,`spatial_corr_3` Functional parcellations with the second and third highest correlations to the spatial maps, and the corresponding correlation values
* `sorted_IC_idx_[atlas].csv` sorted index of the input spatial maps corresponding to the brain networks (artifact-related component indexes are removed)
* `sorted_fnc_[atlas].csv` sorted functional network connectivity (FNC) matrix of the brain networks 
* `nc` folder contains the noisecloud toolbox output. It has the following files:
    * `*.nii` template `nii` files warped into the same space as the input spatial maps.
    * `nc_class_labels.txt` a vector of 0/1 corresponding to the input spatial maps; 0=artifact, 1=network
    * `training/testing_features.csv` contains the training/testing input data features used by the noisecloud toolbox in classification.  

## Result

The following figures are generated using the `./src/example_plot_fnc.m` script. 
You can update the ICA parameter file and autolabeller output folder locations in the above to generate new figures.
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


