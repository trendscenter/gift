**GIFT fMRI Example Data (Resting State)**
### Table of Contents
1. [Introduction](#secIntro)
2. [Data](#secData)
3. [Preprocessing the demo data](#secProcDemo)
4. [Running GIFT](#secGift)
5. [Conclusion](#secConc)
6. [References](#secRef)

# **Introduction** <a name="secIntro"></a>


## Demonstration data for GIFT

Notice Oct 14, 2022: This demo dataset that may be published in November, 2022 does still not exist as it is a work in progress. If you have any questions please email ceierud@gsu.edu.

GIFT is a handy and efficient tool that performs customizable group independent component analysis (GICA) on a study cohort. In this tutorial, we walk you through a typical GIFT analysis explaining the pipeline and its parts. We demonstrate how to use the pipeline on a cohort of 10 males and 10 females from an undisclosed study.


## Different Aspects In Data Processing 

This pipeline processes raw data to the end product, including preprocessing using fmriprep. GICA postprocessing is performed afterwards. This group ICA is guided by Neuromark, a brain atlas derived from a big cohort of fMRI scans, as described in [Du et al. 20201](https://www.sciencedirect.com/science/article/pii/S2213158220302126). At last comes the Dynamic Functional Connectivity step. It calculates connectivity between different brain regions and highlights differences between two groups (healthy controls versus patients).


## Data Not For Research (Disclaimer) 

This data is partitioned to show results with few subjects and is biased and may not be used for research.


# **Data** <a name="secData"></a>


## Raw Data

Raw fMRI and structural T1 data is available if you have time to run demo from scratch. Please contact the authors.

Preprocessed Data

To make the computer processing less daunting to you we have a dataset that has preprocessed all the regular steps ahead.


# **Processing The Demo Data** <a name="secProcDemo"></a>

## `bidsify_neuromark_raw.sh` - Format Dataset According to BIDS

BIDS is a guideline providing a tidy and reproducible way to work with neuroimaging dataset. This script adds a dataset description, participants list and formats the names according to the specification. The BIDS dataset will be used in the downstream applications.

## `run-fmriprep.sh` - Preprocessing the Raw Data 

Data preprocessing includes alignment with the MNI space and different other standartization/artifact removal procedures. We use fMRIprep (Esteban et al., 2019)<sup><a href="#bookmark=id.x3n4vrij65zt">2</a></sup> for data preprocessing. It is possible to run this tool on cluster using e.g. Singularity (or Docker):


```
$ singularity run --cleanenv fmriprep.simg 
    path/to/data/dir path/to/output/dir 
    participant 
    --participant-label label
```


Please refer to [fMRIprep documentation3](https://fmriprep.org/en/1.5.1/index.html) for further information on preprocessing steps and methods. Furthermore, we attach our fMRIprep run script to the present repository under run-fmriprep.sh.

## `smooth_fmriprep_results.sh` and `smooth_subjects.sh` - Smooth the fMRIprep output
For a better ICA reconstruction in GIFT, we smooth the data using a Gaussian kernel. The core functionality is provided by `smooth_subjects.sh`:


```
fslmaths /out/fmriprep/sub-01/func/sub-01_task-mixedgamblestask_run-*1_space-MNI152NLin2009cAsym_desc-preproc_bold.nii.gz  #input file
-kernel gauss 4.2466452  #smoothing kernel details
-fmean #mean type 
/out/fmriprep/sub-01/func/sub-01_task-mixedgamblestask_run-1_space-MNI152NLin2009cAsym_desc-preproc_bold10sm.nii.gz #output path
```
To run `fslmaths`, I start the fmriprep container in `smooth_fmriprep_results.sh` and execute `smooth_subjects.sh` in this container.

This command utilizes fslmaths utility from Freesurfer (available in fmriprep too). Of course, you are welcome to use a different smoothing tool to achieve the same end.

# Running GIFT <a name="secGift"></a>

Now, that the data has been preprocessed, we turn to GIFT to actually carry out the Independent component analysis. GIFT is available in several flavors: as a Docker app, Matlab app with a graphical user interface. In following paragraphs, we look closer at all of the different ways to run GIFT.


## `run-all-gift-neuromark.sh` - Processing Independent Component Analysis Using GIFT-BIDS-App

GIFT-BIDS is available as a Docker container from Docker hub. The flavor of GIFT-BIDS, which we refer to as “regular GIFT”, is available as MATLAB GUI application. 

Please find the source code and Docker link for GIFT-BIDS at [trendscenter/gift-bids (github.com)](https://github.com/trendscenter/gift-bids). You are welcome to report bugs and suggest changes. This repository also includes a small demo on GIFT-BIDS.

GIFT-BIDS does **not **require you to have a MATLAB license.

We launch GIFT-BIDS with the command of the following form: 


```
singularity run --bind <root>/tmp1:/tmp --bind &lt;root>/tmp2:/var/tmp  #bind tmp directories 
--bind <root>/ZN_Neuromark_BIDS:/data  #input data BIDS-formatted directory 
--bind <root>/<output directory>:/output  #output directory 
--bind <root>/cfg:/cfg  #directory with the run configuration, explained below 
<root>/trends_gift-bids.img #singularity GIFT-BIDS container 
/data /output  #pointer to the mounted directories which should be used as in- and output 
participant --participant_label 004 033 111 201  #list of participant IDs to process  
--config /cfg/config_spatial_ica_bids.m  #GIFT run config  
1><log-file> 2>&1  #reroute err and std output to log-file
```


where `<root>` is a path to a root directory. Here, we mount the directories from our system into the Singularity container of GIFT. 

In this command we use Singularity, an alternative virtualization engine to Docker. Docker images and Singularity images are mutually convertible. One can e.g. use docker2singularity for image conversion. Both Docker and Singularity utilize very similar concepts - please consult the respective documentation for more information.

We bind the temporary directories tmp1 and tmp2 as they are required for temporary GIFT files during the run. Furthermore, please note that the input directory _is required to_ have BIDS formatted data. Also, subjects should have IDs of form sub-001, sub-002, sub-003 and so on. This is because GIFT currently outputs subject results in this ID format.

If you have a different form of subject IDs (say, the first is “M87395841”), you might put this ID in participants.tsv stating the new subject_id as “sub-001”. Later on, you could look the old id in this table. You can consult our way of doing this in bidsify_neuromark_raw.sh script in this directory.

Another core part of the analysis is the run configuration.


## GIFT run configuration 

While MATLAB GUI version of GIFT allows to configure the run in the GUI directly, GIFT-BIDS utilizes a pre-written *_config.m file to choose the run mode.

Config file contents are simple “key=value” pairs written in MATLAB. For example, here is an excerpt of config_spatial_ica_bids.m we use in this demo:

```
%% Modality. Options are fMRI and EEG 
modalityType = 'fMRI'; 

%% Data Pre-processing options 
% 1 - Remove mean per time point 
% 2 - Remove mean per voxel 
% 3 - Intensity normalization 
% 4 - Variance normalization
preproc_type = 1; 

%% 'Which ICA Algorithm Do You Want To Use'; 
% see icatb_icaAlgorithm for details or type icatb_icaAlgorithm at the 
% command prompt. 
% Note: Use only one subject and one session for Semi-blind ICA. Also specify atmost two reference function names 
 
% 1 means infomax, 2 means fastICA, etc. 
algoType = 'moo-icar';
```

This signifies the basic information GIFT-BIDS needs to run, as well as choices like ICA algorithm.

To find out more parameter possibilities, check out the configuration file examples under [https://github.com/trendscenter/gift/tree/master/GroupICATv4.0c/icatb/icatb_batch_files](https://github.com/trendscenter/gift/tree/master/GroupICATv4.0c/icatb/icatb_batch_files).


# **Conclusion** <a name="secConc"></a>

We are happy if GIFT-BIDS can deem helpful in your work. Hopefully, this demo made you step closer to utilizing GIFT in your analyses, thus reducing computational burden and processing time.


## Please cite GIFT 

If you have used GIFT in your work, please cite:

 




# References <a name="secRef"></a>

1.  Du, Y., Fu, Z., Sui, J., Gao, S., Xing, Y., Lin, D., ... & Alzheimer's Disease Neuroimaging Initiative. (2020). NeuroMark: An automated and adaptive ICA based pipeline to identify reproducible fMRI markers of brain disorders. _NeuroImage: Clinical_, _28_, 102375.

2.  Esteban, O., Markiewicz, C. J., Blair, R. W., Moodie, C. A., Isik, A. I., Erramuzpe, A., ... & Gorgolewski, K. J. (2019). fMRIPrep: a robust preprocessing pipeline for functional MRI. _Nature methods_, _16_(1), 111-116.

3.  https://fmriprep.org/en/1.5.1/index.html
