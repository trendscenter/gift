# GIFT 
<!-- PLEASE DO NOT EDIT THIS LINE OR LINE BELOW -->
### Group ICA/IVA software (MATLAB) v4.0.5.0
<!-- PLEASE DO NOT EDIT ABOVE THIS LINE -->
![TReNDS](https://trendscenter.org/wp-content/uploads/2019/06/background_eeg_1.jpg)
### Table of Contents
1. [Introduction](#secIntro)
2. [Download](#secDownload)
3. [GIFT BIDS-Apps](#secBids)
4. [Screen Shots](#secScreen)
5. [Toolboxes](#secTools)
	1. [Mancovan](#secToolMan)
	2. [NBiC](#secToolNbic)
6. [Version History](#secVerHist)
---
### Introduction <a name="secIntro"></a>
GIFT is an application supported by the NIH under grant 1RO1EB000840 to Dr. Vince Calhoun and Dr. Tulay Adali. It is a MATLAB toolbox which implements multiple algorithms for independent 
component analysis and blind source separation of group (and single subject) functional magnetic resonance imaging data. GIFT works on MATLAB R2008a and higher. Many ICA algorithms were 
generously contributedby Dr. Andrzej Cichocki. These are also available in Dr. Cichocki's ICALAB toolbox. For any question or comments please contact Vince Calhoun (vcalhoun@gsu.edu) or 
Cyrus Eierud (ceierud@gsu.edu).

Please note that all the toolboxes in GIFT require only MATLAB and not dependent on additional MATLAB toolboxes like Image Processing, Signal Processing, etc. Basic GIFT analysis (without GUI) 
runs on MATLAB R13 and higher. GIFT GUI works on R2008a and higher. 

### Downloads <a name="secDownload"></a>
**GroupICAT**  - Download latest version by clicking the green code button on the upper right on this page and then clone the software using the link and the git clone command in your terminal. Current version of Group ICA. Requires MATLAB R2008a and higher.
#### Stand Alone Versions
[**Windows 64**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/software/stand_alone/GroupICATv4.0c_standalone_Win64.zip) - Compiled on Windows 64 bit OS and MATLAB R2020a. Please see read me text file for more details.\
[**Linux-x86-64**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/software/stand_alone/GroupICATv4.0.3.3_standalone_Linux_x86_64.zip) - Compiled on Linux-x86-64 bit OS and MATLAB R2016b. Please see read me text file for more details.\
[**fMRI Data**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/example_subjects.zip) - Example fMRI datais from a visuomotor paradigm.\
[**Mancovan Sample Data**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/mancova_sample_data.zip) - Sample data to use in mancovan analysis or temporal dfnc analysis.\

[**Complex GIFT**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/software/GroupICATv2.0d_complex.zip) - ICA is applied on complex fMRI data. Please follow the read me text file instructions for doing complex fMRI ICA analysis.\

### GIFT BIDS-Apps <a name="secBids"></a>
If you have your data in BIDS format or you want to run GIFT under a cluster you may want to our GIFT BIDS-Apps [gift-bids](https://github.com/trendscenter/gift-bids). 

### Screen Shots <a name="secScreen"></a>

| ![GIFT](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/images/gift.jpg) |
|:--:|
| Figure 1. Main menu of GIFT|


### Toolboxes <a name="secTools"></a>
#### Mancovan <a name="secToolMan"></a>
Mancovan toolbox is based on the paper (E. Allen, E. Erhardt, E. Damaraju, W. Gruner, J. Segall, R.
Silva, M. Havlicek, S. Rachakonda, J. Fries, R.Kalyanam, A. Michael, J. Turner, T. Eichele, S.
Adelsheim, A. Bryan, J. R. Bustillo, V. P. Clark, S. Feldstein,F. M. Filbey, C. Ford, et al, 2011). This
toolbox works on MATLAB versions greater than R2008a. Features used are subject component
spatial maps, timecourses spectra and FNC correlations. Multivariate tests are done on the features
to determine the significant covariates which are later used in the univariate tests on each feature.
To invoke the toolbox, select â€œMancovanâ€? under â€œToolboxesâ€? menu (Figure 3.2). You could also
invoke toolbox using mancovan_toolbox at the command prompt. Mancovan toolbox (Figure 3.38)
is divided into four parts like create design matrix, setup features, run mancova and display.
#### N-BiC <a name="secToolNbic"></a>
NBiC toolbox is based on the 2020 publication "N-BiC: A Method for Multi-Component and Symptom Biclustering of Structural MRI Data: Application to Schizophrenia" (Md Abdur Rahaman , Jessica A. Turner, Cota Navin Gupta, Srinivas Rachakonda, Jiayu Chen , Jingyu Liu , Theo G. M. van Erp, Steven Potkin, Judith Ford, Daniel Mathalon, Hyo Jong Lee, Wenhao Jiang, Bryon A. Mueller, Ole Andreassen, Ingrid Agartz, Scott R. Sponheim , Andrew R. Mayer, Julia Stephen , Rex E. Jung, Jose Canive, Juan Bustillo, and Vince D. Calhoun). This toolbox works on MATLAB versions greater than R2008a. [Click here for more info](https://github.com/trendscenter/gift/blob/master/GroupICATv4.0c/icatb/toolbox/nbic/README.md).

### Version History<a name="secVerHist"></a>
IcaTbVersion: 4.0.3.5. More information about about the GIFT version history is found at the following link: [GIFT version history](https://github.com/trendscenter/gift/blob/master/doc/updates/README.md) 

