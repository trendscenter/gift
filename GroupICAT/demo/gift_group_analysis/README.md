**GIFT fMRI Example Data (Resting State)**
### Table of Contents
1. [Introduction](#secIntro)
2. [Data](#secData)
3. [Preprocessing the demo data](#secProcDemo)
4. [Running GIFT](#secGift)
5. [Groups differences with MANCOVAN](#secMancovan)
6. [Conclusion](#secConc)
7. [References](#secRef)

# **Introduction** <a name="secIntro"></a>


## Demonstration data for GIFT

Notice Oct 14, 2022: This demo dataset that may be published in November, 2022 does still not exist as it is a work in progress. If you have any questions please email ceierud@gsu.edu.

GIFT is a handy and efficient tool that performs customizable group independent component analysis (GICA) on a study cohort. In this tutorial, we walk you through a typical GIFT analysis explaining the pipeline and its parts. We demonstrate how to use the pipeline on a cohort of 15 males and 15 females from an undisclosed study.


## Different Aspects In Data Processing 

This pipeline processes raw data to the end product, including preprocessing using fmriprep. GICA postprocessing is performed afterwards. This group ICA is guided by Neuromark, a brain atlas derived from a big cohort of fMRI scans, as described in [Du et al. (2020)](#refDu2020). Finally one may postprocess the ICA using MANCOVAN or Dynamic Functional Connectivity step. It calculates connectivity between different brain regions and highlights differences between two groups (healthy controls versus patients).


## Data Not For Research (Disclaimer) 

This data is partitioned to show results with few subjects and is biased and may not be used for research.


# **Data** <a name="secData"></a>


## Raw Data

Raw fMRI and structural T1 data is available if you have time to run demo from scratch using a different dataset that does not support the MANCOVAN features, but it supports GIFT blind ICA and other non- MANCOVAN ICA for up to 10 subjects. Download the ds005.tar dataset from https://drive.google.com/drive/folders/0B2JWN60ZLkgkMGlUY3B4MXZIZW8?resourcekey=0-EYVSOlRbxeFKO8NpjWWM3w into the demo directory. 

If you want to step to the "Running GIFT" section you may start with preprocessed data found [here](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/demo_input3neuromark.zip).

If you just want to run the last step you may choose "root/demo_input3neuromark/gift_out/neuromark_ica_parameter_info.mat" as the parameter file, following instructions under section "Groups differences with MANCOVAN".

# **Processing The Demo Data** <a name="secProcDemo"></a>

## `bidsify_neuromark_raw.sh` - Format Dataset According to BIDS

Brain imaging data structure (BIDS) is a guideline providing a reproducible way to work with neuroimaging dataset. This script adds a dataset description, participants list and formats the names according to the specification. The BIDS dataset will be used in the downstream applications. This script is optional to run the GIFT demo as we also provide the upstream dataset after the time consuming fMRIprep processing.

## `run-fmriprep.sh` - Preprocessing the Raw Data 

Data preprocessing includes normalization from native space into the MNI space and different other standartization/artifact removal procedures. We use fMRIprep ([Esteban et al., 2019](#refEsteban2019)) for data preprocessing. Also this step is optional (as above) as you may save time skipping to the "Running GIFT" section (using the dataset that already has completed the time consuming fMRIprep processing). Please refer to [fMRIprep documentation](https://fmriprep.org/en/1.5.1/index.html) for installation and detailed preprocessing steps and methods. fMRIprep may be run using a container (Singularity or Docker) and this example shows how to run Singularity:


```
$ singularity run --cleanenv fmriprep.simg 
    path/to/data/dir path/to/output/dir 
    participant 
    --participant-label label
```


Furthermore, we attach our fMRIprep run script to the present repository under run-fmriprep.sh.

## Smooth the fMRIprep output
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
If you did not run the previous steps you may download the dataset, containing both input and output files for GIFT. Then you may run GIFT with your own output directory. The input and output files are found [here](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/demo_input3neuromark.zip).<br>
Only the preprocessed fMRI is found [here](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/demo_input3neuromark_onlysubj.zip).

Now, that the data has been preprocessed, we turn to GIFT to actually carry out the Independent component analysis. GIFT is available in several flavors: as a Docker app, Matlab app with a graphical user interface. In following paragraphs, we look closer at different ways to run GIFT.


## Processing Independent Component Analysis Using GIFT-BIDS-App

GIFT-BIDS is available as a Docker container from Docker hub. The flavor of GIFT-BIDS, which we refer to as “regular GIFT”, is available as MATLAB GUI application. 

Please find the source code and Docker link for GIFT-BIDS at [trendscenter/gift-bids (github.com)](https://github.com/trendscenter/gift-bids). You are welcome to report bugs and suggest changes. This repository also includes a small demo on GIFT-BIDS.

GIFT-BIDS does **not **require you to have a MATLAB license.

We launch GIFT-BIDS with the command of the following form (run-all-gift-neuromark.sh): 


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

If you have a different form of subject IDs (say, the first is “M87395841”), you put this ID in participants.tsv stating the new subject_id as “sub-XXX”. Later on, you could look the old id in this table. You can consult our way of doing this in bidsify_neuromark_raw.sh script in this directory. If you don't do so, you could have a hard time matching subject IDs as output by GIFT with your original ones later in the analysis.

Another core part of the analysis is the run configuration.

*GIFT run configuration*
Too highlight the different possibilities the gift-bids app has, showing other possibilities a small background of the different batch files flavors exists. While MATLAB GUI version of GIFT allows to configure the run in the GUI directly, GIFT-BIDS utilizes a pre-written *_config.m file to specify all the necessary properties (key values) before the run.

The config file contains “key=value” pairs written in MATLAB. For example, here is an excerpt of config_spatial_ica_bids.m we use in this demo:

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

# Groups differences with MANCOVAN <a name="secMancovan"></a>

To finally find the differences between the groups visualize them, we will run the MANCOVAN toolbox, included in GIFT. MANCOVAN allows to compare individual ICs and subjects using statistical tests. In this example we will use the GIFT GUI, where you need to start MATLAB and then in the MATLAB command line enter a couple of lines to start gift:
-  addpath(genpath('/my/gift/path'));
-  gift

To launch MANCOVAN, select it from the main GIFT toolbox menu (Fig. 1). ![images/image1.png](images/image1.png)<br> 
Figure 1.

The main MANCOVAN window appears. 
We need to create design matrix first. Click on the corresponding box (Fig. 2). The window will prompt you for a parameter file.
It is located in the output directory of the GIFT analysis (output from gift-bids app). ![images/image2.png](images/image2.png)Figure 2.

In the next step, traverse down to a MANCOVAN output directory (that you perhaps have previously made) and press "." in the right pane to chose the directory and press "OK" in the bottom to confirm.

Next, the test configuration window appears. We will be using a "2-sample t-test" with no covariates (Fig. 3). <a name="image3"></a>
![images/image3.png](images/image3.png)Figure 3.

In the screen that appears, define groups. Give a name to a group (e.g., males), and hold "CTRL" to select multiple subjects of the male group (e.g.: 1, 3, 4, 5, 10, 12, 14, 15, 17, 18, 21, 23, 25, 29, 30, please note that Fig. 4 is slightly wrong). Click "OK" when all the subjects of a group have been highlighted. 
![images/image4.png](images/image4.png)Figure 4.

Now, define the second, female groups, analogously as the first group, but with different subjects (e.g.: 2, 6, 7, 8, 9, 11, 13, 16, 19, 20, 22, 24, 26, 27, 28). ![images/image5.png](images/image5.png)Figure 5.

Click "OK". Return to the "Setup MANCOVAN Design" window ([Fig. 3](#image3)). Click "Create..." in the bottom of the window.
It gets us back to the main menu (Fig. 6). The design matrix is set up. It is time to set up features. 
Click the [Setup Features] button. ![images/image6.png](images/image6.png)Figure 6.

A prompt to select mancovan parameter file appears (Fig. 7). It is located in the output directory 
you have appointed in the previous steps. ![images/image7.png](images/image7.png)Figure 7.

After selecting the parameter file, MANCOVAN Setup Analysis appears. 
We will be using FNC correlation without lags (Fig. 8). Then press the [+] button to add components.
![images/image8.png](images/image8.png)Figure 8.

In this demo we will use the Neuromark template (Figures 9-10).
![images/image9.png](images/image9.png)<br> 
Figure 9. 

![images/image10.png](images/image10.png)<br> 
Figure 10.

In accordance with Fig. 11, name the main network and select its subcomponents, selecting multiple indices (hold CTRL to select multiple manually, or SHIFT to select a range). To see each subnetwork you have selected you may press [S...] at the bottom. Which subnetwork that belongs to the main network is found in Fig. 9. 
![images/image11.png](images/image11.png)Figure 11.

Press "Done" to confirm the network name, having its subnetworks. The defined network is now enlisted. Repeat for all the components you wish to define.
![images/image12.png](images/image12.png)Figure 12.

We are free to change the P-Value threshold and TR accordingly. To define Number of components for each vector manually, proceed as following:
first, tick "Autoselect No. of components..." ![images/image13.png](images/image13.png)Figure 13.

Then, untick "Autoselect No. of components...". An input frame appears (Fig. 14). Enter 15 as the default number of components (2, as shown in Fig. 14 may also work). ![images/image14.png](images/image14.png)Figure 14.

Click "Run" (Fig. 14) in the bottom of the window. MANCOVAN loads subjects. ![images/image15.png](images/image15.png)Figure 15.

When loading is complete, we are ready to run MANCOVAN (Fig. 16). ![images/image16.png](images/image16.png)Figure 16.<br> 

Press "Run MANCOVAN" in the main menu (Fig. 16). We do not want to remove nuisance covariates (Fig. 17). ![images/image17.png](images/image17.png)Figure 17.<br> 

Select MANCOVAN parameter file as done previously (in the projected MANCOVAN output directory). MANCOVAN prints the output
to the MATLAB console. ![images/image18.png](images/image18.png)Figure 18.

After MATLAB ceases to display "Processing..." marker in the left bottom part, we can display the results (Fig. 19). ![images/image19.png](images/image19.png)Figure 19.

Display univariate results (Fig. 20). ![images/image20.png](images/image20.png)Figure 20.<br> 

Define T-Threshold to be 1.0, positive and negative image values. We will not do any multiple comparisons correction (too few subjects), leave fALFF defaults and display the connectogram (Fig. 21). ![images/image21.png](images/image21.png)Figure 21.<br> 

Connectogram appears. We see significant differences in the functional connectivities between components (Fig. 22). ![images/image22.png](images/image22.png)Figure 22.

Another figure gives a hint on how significantly different components pairs differ between males and females in a matrix (Fig. 23).
![images/image23.png](images/image23.png)Figure 23.

The final figure (Fig. 24) the difference between groups for each main network. It looks 
like in our setup, Default Mode network and Cognitive-Control network might have significantly different functional 
connectivity patterns between males and females. ![images/image24.png](images/image24.png)Figure 24.


# **Conclusion** <a name="secConc"></a>

We are happy if GIFT-BIDS and adjacent toolboxes can deem helpful in your work. Hopefully, this demo helped you a step closer to utilizing GIFT in your analyses, thus reducing computational burden and processing time.


# References <a name="secRef"></a>

1.  Du, Y., Fu, Z., Sui, J., Gao, S., Xing, Y., Lin, D., ... & Alzheimer's Disease Neuroimaging Initiative. (2020). NeuroMark: An automated and adaptive ICA based pipeline to identify reproducible fMRI markers of brain disorders. _NeuroImage: Clinical_, _28_, 102375. <a name="refDu2020"></a> [Click here for article](https://www.sciencedirect.com/science/article/pii/S2213158220302126)

2.  Esteban, O., Markiewicz, C. J., Blair, R. W., Moodie, C. A., Isik, A. I., Erramuzpe, A., ... & Gorgolewski, K. J. (2019). fMRIPrep: a robust preprocessing pipeline for functional MRI. _Nature methods_, _16_(1), 111-116. <a name="refEsteban2019"></a>
