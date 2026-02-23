# GIFT 
<!-- PLEASE DO NOT EDIT THIS LINE OR LINE BELOW -->
### Group ICA/IVA software (MATLAB) v4.0.6.32
<!-- PLEASE DO NOT EDIT ABOVE THIS LINE -->
![TReNDS](https://trendscenter.org/wp-content/uploads/2019/06/background_eeg_1.jpg)
### Announcements
- Apple Silicon! If your Apple Silicon computer (M1, M2, M3 or M4) gives you errors you may need to unquarantine your mex files by doing following from a Mac terminal:
cd /my/software/folder/gift/GroupICAT/icatb && find . -iname "\*.mexmaca64" -exec xattr -d com.apple.quarantine {} \\; && find . -iname "\*.mexmaci64" -exec xattr -d com.apple.quarantine {} \\;
- GIFT v4.0.5.14 (10/31/2024) slightly modified the dFNC and MANCOVAN processing order, now having the following steps: 1) detrending, 2) regressing out confounds, 3) despiking, 4) filtering.
### Table of Contents
1. [Introduction](#secIntro)
2. [Download](#secDownload)
3. [GIFT BIDS-Apps](#secBids)
4. [Screen Shots](#secScreen)
5. [Version Compatability](#verComp)
6. [Toolboxes and Features added to GIFT](#secTools)
	1. [Mancovan](#secToolMan)
 	2. [EEGIFT](#secToolEEGIFT)
	3. [NBiC](#secToolNbic)
 	4. [Noise Cloud](#secToolNoise)
  	5. [Autolabeller](#secToolAutolabeller) 
7. [Documentation/Manual](#manual)
8. [FAQ](#faq)
9. [Version History](#secVerHist)
10. [Publications](#pubs)
---
### Introduction <a name="secIntro"></a>
GIFT is an application, originally supported by NIH grant 1RO1 EB000840 to Dr. Vince Calhoun and Dr. Tulay Adali and has been continuous supported by NIH and NSF. The MATLAB application implements multiple algorithms for independent component analysis and blind source separation of group (and single subject) functional magnetic resonance imaging data. GIFT has been downloaded 18581 times (as of 7/6/24) by researchers world wide. For question or comments please contact Vince Calhoun (vcalhoun@gsu.edu) or Cyrus Eierud (ceierud@gsu.edu).

### Downloads <a name="secDownload"></a>
**GroupICAT**  - Download latest version by clicking the green code button on the upper right on this page and then clone the software using the link and the git clone command in your terminal. Current version of Group ICA. Requires MATLAB R2008a and higher.
#### Stand Alone Versions
[**Windows 64**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/software/stand_alone/GroupICATv4.0c_standalone_Win64.zip) - Compiled on Windows 64 bit OS and MATLAB R2020a. Please see read me text file for more details.\
[**Linux-x86-64**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/software/stand_alone/GroupICATv4.0.3.3_standalone_Linux_x86_64.zip) - Compiled on Linux-x86-64 bit OS and MATLAB R2016b. Please see read me text file for more details.\
[**fMRI Data**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/example_subjects.zip) - Example fMRI datais from a visuomotor paradigm.
[**Mancovan Sample Data**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/data/mancova_sample_data.zip) - Sample data to use in mancovan analysis or temporal dfnc analysis.

[**Complex GIFT**](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/gift/software/GroupICATv2.0d_complex.zip) - ICA is applied on complex fMRI data. Please follow the read me text file instructions for doing complex fMRI ICA analysis.

### GIFT BIDS-Apps <a name="secBids"></a>
If you have your data in BIDS format or you want to run GIFT under a cluster you may want to our GIFT BIDS-Apps [gift-bids](https://github.com/trendscenter/gift-bids). 

### Screen Shots <a name="secScreen"></a>

| ![GIFT](https://github.com/trendscenter/gift/blob/master/doc/web/img/20240705Gift4Ims.png) |
|:--:|
| Figure 1. GIFT Sceenshots|

### Version Compatability <a name="verComp"></a>
All the toolboxes in GIFT require only MATLAB and not dependent on additional MATLAB toolboxes like Image Processing, Signal Processing, etc. Basic GIFT analysis (without GUI) runs on MATLAB R13 and higher. GIFT GUI works on R2008a and higher. Please see below for specific details related to toolboxes in GIFT:

- MANCOVAN runs on MATLAB R2008a and higher. From R2012b onwards, Optimization toolbox is not required to compute t-threshold based on distribution of voxelwise t-stats. There is an option to use Z-threshold or select mask if Optimization toolbox is not installed on MATLAB versions less than R2012b.
- Temporal and Spatial dFNC runs on MATLAB R2008a and higher.


### Toolboxes and Features added to GIFT <a name="secTools"></a>
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
#### EEGIFT <a name="secToolEEGIFT"></a>
[Click here for EEGIFT documentation](https://trends-public-website-fileshare.s3.amazonaws.com/public_website_files/software/eegift/docs/v1.0c_EEGIFT_Walk_Through.pdf)

#### N-BiC <a name="secToolNbic"></a>
NBiC toolbox is based on the 2020 publication "N-BiC: A Method for Multi-Component and Symptom Biclustering of Structural MRI Data: Application to Schizophrenia" (Md Abdur Rahaman , Jessica A. Turner, Cota Navin Gupta, Srinivas Rachakonda, Jiayu Chen , Jingyu Liu , Theo G. M. van Erp, Steven Potkin, Judith Ford, Daniel Mathalon, Hyo Jong Lee, Wenhao Jiang, Bryon A. Mueller, Ole Andreassen, Ingrid Agartz, Scott R. Sponheim , Andrew R. Mayer, Julia Stephen , Rex E. Jung, Jose Canive, Juan Bustillo, and Vince D. Calhoun). This toolbox works on MATLAB versions greater than R2008a. [Click here for more info](https://github.com/trendscenter/gift/blob/master/GroupICAT/icatb/toolbox/nbic/README.md).
#### Noise Cloud <a name="secToolNoise"></a>
Noise cloud uses both spatial and temporal features to identify noise/artifact components from the specified components. More information about Noisecloud will soon come here...
#### Autolabeller <a name="secToolAutolabeller"></a>
Autolabeller info will soon come here ...

### Documentation/Manual<a name="manual"></a>
[Click here for link to manual in PDF format](https://github.com/trendscenter/gift/blob/master/doc/gica_manual.pdf) <br>
[Click here for link to manual in Word format](https://github.com/trendscenter/gift/blob/master/doc/gica_manual.docx)

### FAQ<a name="faq"></a>
[Click here for FAQ](https://github.com/trendscenter/gift/blob/master/doc/web/faq/README.md)

### Version History<a name="secVerHist"></a>
More information about about the GIFT version history is found at the following link: [GIFT version history](https://github.com/trendscenter/gift/blob/master/doc/web/updates) 

### Publications<a name="pubs"></a>
[Click here for a list of GIFT related publications](https://github.com/trendscenter/gift/blob/master/doc/web/publications/README.md) 

