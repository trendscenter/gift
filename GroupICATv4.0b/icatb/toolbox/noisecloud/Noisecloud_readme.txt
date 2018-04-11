Noise cloud tool is used to detect noise components. It is assumed that SPM is installed on MATLAB path. Inputs are training set (known subject component spatial maps and timecourses) and 
testing set (unknown subjects that needs to be classified). You need to select labels associated with the training set. To use the tool, select GIFT -> Toolboxes -> Noisecloud:

1. After you selected the output directory for the analysis, a GUI will open which contains two panels like testing and training.

2. Select the training components, you could select multiple subject component nifti files (*sub*comp*nii). Number of total components (subjects x components) must match when selecting 
training timecourses (*sub*time*nii). Enter TR of the experiment in seconds followed by a labels text file. Labels text file must be ascii and must contain binary numbers (1 or 0) where 1 
corresponds to noise and 0 for component network. Length of labels vector in the file must match the number of total components selected. For example, if you selected 10 subjects each 
containing 50 components, labels file must contain a vector of length 1 x 500.

3. Enter the testing data (component maps and timecourses) which needs to be classified. At the end of the analysis, confusion matrix is shown in the graphical window. 
Also results are saved in file noise_cloud_results.mat. Variables description are as follows:
	a. class_labels variable in the MAT file which contain the flags associated with the testing components (Noise or network)
	b. result_nc_classifier - Classifier built using the training data
	c. fit_mdl - Model fit
	d. training_opts - Training options
	e. testing_opts - Testing options


Alternatively, you could invoke noise cloud from the command line using the function noisecloud_run. Please see help noisecloud_run or noisecloud_demo.m for more details.