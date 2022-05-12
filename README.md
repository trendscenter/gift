# GIFT 
###Group ICA/IVA software (MATLAB)
![TReNDS](https://trendscenter.org/wp-content/uploads/2019/06/background_eeg_1.jpg)
### Table of Contents
1. [Introduction](#secIntro)
2. [Download](#secDownload)
3. [Screen Shots](#secScreen)
4. [Toolboxes](#secTools)
a. [Autolabeller](#secPluAuto)
5. [Version History](#secVerHist)
---
### Introduction <a name="secIntro"></a>
GIFT is an application supported by the NIH under grant 1RO1EB000840 to Dr. Vince Calhoun and Dr. Tulay Adali. It is a MATLAB toolbox which implements multiple algorithms for independent 
component analysis and blind source separation of group (and single subject) functional magnetic resonance imaging data. GIFT works on MATLAB R2008a and higher. Many ICA algorithms were 
generously contributedby Dr. Andrzej Cichocki. These are also available in Dr. Cichocki's ICALAB toolbox. For any question or comments please contact Vince Calhoun (vcalhoun@gsu.edu) or 
Cyrus Eierud (ceierud@gsu.edu).

Please note that all the toolboxes in GIFT require only MATLAB and not dependent on additional MATLAB toolboxes like Image Processing, Signal Processing, etc. Basic GIFT analysis (without GUI) 
runs on MATLAB R13 and higher. GIFT GUI works on R2008a and higher. 

### Downloads <a name="secDownload"></a>
**GroupICAT v4.0c**  - Download by clicking the green clone button on the upper left on this page and then clone the software using hte link using the git clone command in your terminal. Current version of Group ICA. Requires MATLAB R2008a and higher.
    **Stand Alone Version**
        [**Windows 64**](https://trendscenter.org/trends/software/gift/software/stand_alone/GroupICATv4.0c_standalone_Win64.zip) - Compiled on Windows 64 bit OS and MATLAB R2020a. Please see read me text file for more details.
        [**Linux-x86-64**](https://trendscenter.org/trends/software/gift/software/stand_alone/GroupICATv4.0c_standalone_Linux_x86_64.zip) - Compiled on Linux-x86-64 bit OS and MATLAB R2016b. Please see read me text file for more details.
    [**fMRI Data**](https://trendscenter.org/trends/software/gift/data/example_subjects.zip) - Example fMRI datais from a visuomotor paradigm.
    [**Mancovan Sample Data**](https://trendscenter.org/trends/software/gift/data/mancova_sample_data.zip) - Sample data to use in mancovan analysis or temporal dfnc analysis.

[**Complex GIFT**](https://trendscenter.org/trends/software/gift/software/GroupICATv2.0d_complex.zip) - ICA is applied on complex fMRI data. Please follow the read me text file instructions for doing complex fMRI ICA analysis.


### Screen Shots <a name="secScreen"></a>

| ![GIFT](https://trendscenter.org/trends/software/gift/images/gift.jpg) |
|:--:|
| Figure 1. Main menu of GIFT|


### Toolboxes <a name="secTools"></a>
#### Mancovan <a name="secPlugAuto"></a>
Mancovan toolbox is based on the paper (E. Allen, E. Erhardt, E. Damaraju, W. Gruner, J. Segall, R.
Silva, M. Havlicek, S. Rachakonda, J. Fries, R.Kalyanam, A. Michael, J. Turner, T. Eichele, S.
Adelsheim, A. Bryan, J. R. Bustillo, V. P. Clark, S. Feldstein,F. M. Filbey, C. Ford, et al, 2011). This
toolbox works on MATLAB versions greater than R2008a. Features used are subject component
spatial maps, timecourses spectra and FNC correlations. Multivariate tests are done on the features
to determine the significant covariates which are later used in the univariate tests on each feature.
To invoke the toolbox, select “Mancovan” under “Toolboxes” menu (Figure 3.2). You could also
invoke toolbox using mancovan_toolbox at the command prompt. Mancovan toolbox (Figure 3.38)
is divided into four parts like create design matrix, setup features, run mancova and display.

### Version History<a name="secVerHist"></a>
Click the following link for the GIFT version history: [GIFT version history](https://trendscenter.org/trends/software/gift/version_history.html) 

