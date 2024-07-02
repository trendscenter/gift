# Frequently Asked Questions
Some of the frequently asked questions are givenbelow:<br>
[Question 1: What preprocessing steps should be done before analyzing a single subject/ single session in GIFT?](#q1)<br>
[Question 2: What is assumed about the image intensity scale? Do images need to be scaled in a certain way, or conform to a certain format?](#q2)<br>
[Question 3: How do I change the defaults?](#q3)<br>
[Question 4: How do I handle out of memory errors?](#q4)<br>
[Question 5: What do I do if I see the following error dialog boxes: "functional data not selected", "out of memory error for estimating components", and "not a valid number for PC"?](#q5)<br>
[Question 6: Are there any guidelines for choosing the number of components?](#q6)<br>
[Question 7: Are there any other reasons to change the default settings for Infomax (other than the "Select Extended 0/1" option)?](#q7)<br>
[Question 8: What is done during the "Group ICA" step for a single subject single session analysis?](#q8)<br>
[Question 9: Are there any benefits to running subjects in a group vs. running them individually?](#q9)<br>
[Question 10: A single subject might be reduced from 53*63*34*220 to 53*63*34*50." During the data reduction step, are specific time points taken out or averaged?](#q10)<br>
[Question 11: What happens if the calculate ICA step doesn't converge during the analysis?](#q11)<br>
[Question 12: How does the SPM2/SPM5/SPM8 design matrix affect ICA?](#q12)<br>
[Question 13: Which components will be displayed during sorting when component explorer visualization method is used?](#q13)<br>
[Question 14: How do I view all subject time courses at once?](#q14)<br>
[Question 15: How to do I quickly do analyses using different parameters like changing the mask or selecting different numbers for data-reduction steps with the same functional data.](#q15)<br>

---

### Question 1: What preprocessing steps should be done before analyzing a single subject/ single session in GIFT?<a name="q1"></a>

Solution: It depends upon what you are doing. Typically the approach used is to preprocess your data the same way you would for, say, an fMRI analysis of another type (e.g. SPM). In this case, the functional data are spatially realigned,  normalized (if desired) and smoothed in order to run ICA. However, if your goal is to examine potential artifacts in your data, then you may not want to smooth, since the artifacts can be "hidden" by the smoothing process.

### Question 2: What is assumed about the image intensity scale? Do images need to be scaled in a certain way, or conform to a certain format?<a name="q2"></a>

Solution: Image intensity need not be scaled and images should be in 3D/4D analyze or 3D/4D Nifti format. All the images should have the same voxel size but the time points or number of images can be different for subjects.

### Question 3: How do I change the defaults?<a name="q3"></a>

Solution: The defaults are specified in icatb_defaults.m file. Please see the appendix section of the GIFT walk-through or manual. You could also change defaults through GUI using groupica defaults at the MATLAB prompt.

### Question 4: How do I handle out of memory errors?<a name="q4"></a>

Solution: Out of memory errors are common for multivariate approaches like ICA. Thus you want to minimize, wherever possible, the amount of needed memory. For example, don't reslice your data into voxels smaller than the acquired voxel size (since this doesn't add information and just increases memory). If you have very high resolution data, you may even in some extreme cases need to resample your data or work on a subset of the data. Finally, you will benefit greatly by closing all existing MATLAB applications and running MATLAB with out its Java Virtual Machine (JVM) . MATLAB with out JVM can be run by entering command matlab -nojvm at the DOS or UNIX prompt (or change your windows shortcut to include the -nojvm flag so it reads, e.g. "..\matlab.exe -nojvm". This will reduce the amount of memory MATLAB needs to run.

### Question 5:What do I do if I see the following error dialog boxes: "functional data not selected", "out of memory error for estimating components", and "not a valid number for PC"?<a name="q5"></a>

Solution: The ICA process is not interrupted with these errors. The explanation for each error is given below:

What do I do when the "functional data not selected" error dialog box appears?
Select the functional data by clicking on push button Select and selected options will be displayed with options 'Yes' and 'No'.
What do I do when the "out of memory" error dialog box appears while estimating components?
This doesn't terminate the setup ICA process, it just means there is not enough memory available to estimate the number of components. You can continue entering parameters and press Done button when finished.

### Question 6: Are there any guidelines for choosing the number of components?<a name="q6"></a>

Solution: The number of components can be estimated by selecting the estimate dimensionality step in setup ICA. If the number of components are estimated to be very high (like 100 components for 200 images) then around 20 to 30 components should give a reasonable answer.

### Question 7: Are there any other reasons to change the default settings for Infomax (other than the "Select Extended 0/1" option)?<a name="q7"></a>

Solution: Probably not, unless the algorithm fails to converge or you are interested in algorithmic comparisons.

### Question 8: What is done during the "Group ICA" step for a single subject single session analysis?<a name="q8"></a>

Solution: For single subject single session analysis, there is only a single data reduction stage, and you cannot compute the group stats (e.g. T-maps).

### Question 9: Are there any benefits to running subjects in a group vs. running them individually?<a name="q9"></a>

Solution: It is possible to run them individually, however it can be difficult to know how to compare components (since they do not come out in any particular order and there may be different components revealed in different subject). It is not as straightforward to do group analyses with ICA as it is for, say, a model-based approach. For more details please refer to the 2001 Group ICA paper in Human Brain Mapping. A related question is whether to run one or two ICA analyses on different groups. It is recommended to run one large ICA and compare the loading parameters after the analysis, or do group ICA on the controls followed by spatially constrained ICA or spatiotemporal regression to compute output for the rest of the data. In summary, the advantages of running subjects in a group over running them individually are as follows:

Statistics are performed for a component over subjects and sessions.
A selected component can be visualized for different subjects and sessions using Subject Explorer visualization method in display GUI.
Components can be compared for different subjects and sessions based on temporal sorting by selecting all data sets. You can visualize which subject session component is best correlated with the regressors selected by clicking the time course window of the component.

### Question 10: A single subject might be reduced from, say, 53*63*34*220 to 53*63*34*50." During the data reduction step, are specific time points taken out or averaged?<a name="q10"></a>

Solution: During the data-reduction step time points are not averaged but calculated using Principal Component Analysis (PCA). PCA finds the orthogonal components and the number of orthogonal components corresponds to the number 50 in the example. This typically corresponds to more than 99% of the variance in the data.

### Question 11: What happens if the calculate ICA step doesn't converge during the analysis?<a name="q11"></a>

Solution: You can select the calculate ICA step instead of running all the group ICA steps. After the ICA step converges, you have to select back reconstruction, calibrate components and group stats steps.

### Question 12: How does the SPM2/SPM5/SPM8 design matrix affect ICA?<a name="q12"></a>

Solution: SPM design matrix doesn't effect ICA process or run analysis step in GIFT (unless you are using the semi-blind ICA method). The SPM design matrix is used to sort components temporally. A detailed explanation of how SPM design matrix is used is given in appendix section of the Walk-through or in the HTML help manual (When you click Help button in the toolbox, the HTML help manual will open).

### Question 13: Which components will be displayed during sorting when component explorer visualization method is used?<a name="q13"></a>

Solution: Components displayed depends on the type of sorting used.
 
Temporal sorting
Time courses displayed are the time courses selected during sorting. Spatial maps are the maps selected in the viewing set of the component explorer.
Spatial Sorting
Spatial maps and time courses plotted are the components selected for sorting.

### Question 14: How do I view all subject time courses at once?<a name="q14"></a>

Solution: All subject time courses can viewed at once by selecting kurtosis as sorting criteria and selecting all data sets in component set to sort. You can view all the time courses by clicking on SplitTimecourses button on the expanded view of the time course figure.

### Question 15: How do I quickly do analyses using different parameters like changing the mask or selecting different numbers for data-reduction steps with the same functional data?<a name="q15"></a>

Solution: You can do different analyses with the same functional data by copying the parameter file and the Subject file to a different directory. After copying the files, use the setup ICA GUI to change the parameters for the new analysis by specifying the same prefix for the output files you have selected before (If the parameter file name is test_ica_parameter.mat then test is the prefix). You need not to select the subjects again.
