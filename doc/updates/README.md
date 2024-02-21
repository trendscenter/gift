# **GIFT Updates**

## GIFT 4.0.5.0 (Feb 20, 2024):
- Changed name of Neuromark template to Neuromark_fMRI_1.0.nii
- Added Neuromark fMRI 2.0
    - Multiorder constrained maps supplied
- Implemented Datavis
- GIFT version displayed in main window
- Added automatic slicing feature for both fMRI and SBM
- Implemented SSB SWPC option by Ashkan
- Default&ICV mask optoin to remove eyeballs
- Implemented NBIC toolbox, including import function and example script
- GIFT MATLAB Runtime Compiler
    - Reports working
    - Removed rspm_progress_bar in deployed cases (preventing a crash)
- Modified report for bidsapp
- MANCOVAN
    - F-stat option when using univariate results (icatb_defaults)
	- contrasts label excludes nan in univariate results
	- ANOVA fixes
- New demo data for GIFT (Stefan Dvoretski )
- Fixed report generation for spatial constraint ica
- Fix of "out of range" error when inconsistency between comp_network_names, sesInfo.numComp and number of subjects appears
- Batch fix for head motion variables
- GIFT will pick the batch file TR for DFNC processing if if TR is different in the parameter file
- Corrected icasso figures plots out as they should when clicking button [Results Summary]
- T-test added to SBM
- ROI based FC stats summary is now saved
- Plotting components faster
- Trilinear interpolation for the display of function on the structural template (icatb_default.m)
- Implemented cEBM to GIFT
- Added colorbar to icatb_overlayImages
- For batch fix so folder is created even if extra slash is supplied
- Fix for matlab2023b bug
- Default&icv mask for GIFT when run on server installation
- Added average mask for fMRI
- Fixed report generation for Windows platform
- ROI-voxel option in dFC tool fixed so voxel maps are not prevented
- groupica.m now handles multiple arguments when called in batch

 
***

## Miscellaneous enhancements & fixes to GroupICAT v4.0c (Feb 24, 2022):
- Using double precision to avoid any errors in the spatial chronnectome when using default despike option
- Display results structure is added in nipype model file, including network summary options
- Gig-ica algorithm name is changed to MOO-ICAR
- Subject ICA loadings are generated when using algorithms like Infomax, fast ICA, etc
- Default template is changed to ch2bet_3x3x3mm.nii
- Evaluation criteria for estimating clusters like daviesbouldin and ray turi
- icatb_nan_mT function handles multiple contrasts
- Kurtosis graphs y limits are changed separately for the spatial maps and timecourses
- Subject loading coefficients are written
- Updated mex binaries, including for SPM12 
- Option added, to remove components from testing data-sets in the noise cloud toolbox
- Timecourses entered as row vectors are internally converted to column vectors
- Merge analysis only uses timecourse information if spatial components information is not present
- Added INTERP_VAL interpolation when resizing images (icatb_defaults.m)
- GIFT now creates output directory when saving concatenated component timecourses
- Nans are used in timecourses or fnc correlations if subject back-reconstucted files are missing
- Added modified getStateCorrs sub-function (originally found in icatb_post_process_dfnc.m) to roi-based dFC post-processing. Adds call to getStateCorrrs after kmeans clustering is run on the full input dataset (all subjects and windows). This provides users with two additional fields (corrs_state and states) in clusterInfo which are needed for group or individual subject analyses.
- Added options to turn off mutual information and kurtosis in reports
- MOO-ICAR option is changed to use reference file names instead of reference data to handle large number of references
- icatb_save ica data is handled to use single component
- Post process timecourses is saved incrementally
- Added options to compute aggregate spectra, fnc, etc in post process step to speedup display results. Files are saved individual subject-wise instead of one big file.
- Options are provided to compute aggregate fnc, spectra, etc in post-process step. Individual subject files are saved instead of one big file.
- "chkSize undefined value", "too many input args", "Brace indexing is not supported for this variable type", missing field "postprocess" and "Error using reshape" errors fixed
- Anisotropic template is resliced to isotropic
- Tall array DFNC option added when using kmeans
- Ratio to interpolate is computed once and used in resampling timecourses when the TRs are different across subjects
- r_to_z function is used instead of atanh
- Added IVA-L-SOS-Adaptive algorithm
    - iva second order is updated to use only weight change as the stopping criteria
	- Options to initialize weights using IVA-G
- Colormap is used in figure property instead of calling colormap function
- Option to store FNC matrices displayed in field wfcInfo.display_info.FNC
- Option to replace gig-ica with MOO-ICAR to read batch file inputs from previous version
- Options to initial centroids as user input in standard dfnc
- Supporting coregistering files when the format is nifti gzip
- Report generator allows results structure in mat file  and is opened in background mode when using GUI to handle empty plots
- options are added in defaults to write stats info (mancova) and spectra. 
    - Spm stats, calculate stats and single trial amplitude are removed from options dropdown box
- Warning message related to eigs function is fixed
- Decentralized option in mancova
- Option is provided to use spectra options like Npoint FFT and bins when TIMECOURSE_POSTPROCESS.spectra.option is set to 2
- Added decentralized mancova options
- Input Kmeans centroids are back-projected on to the data to find centroids on the new data-set
- Power spectra is computed using pwelch
- Fixed gii file issue which gives error no private field
- Added algorithm table
- Use gray instead of white color for separating cells in the matrix plot
- FNC network plot is used when network names are passed in the univariate results
 
***

## GroupICAT v4.0c (Oct 10, 2020):
- We upgraded some tools like adding neuromark template in constrained ica, display summary tools for source based morphometry and mancovan toolbox in the stand alone version of gift. Docker for group ica is now available. Please download tools at https://trendscenter.org/software/gift/. Docker can also be accessed using https://github.com/trendscenter/gift-bids.
- New GUI is provided to run automated ICA algorithms like MOO-ICAR (previously GIG-ICA) and Constrained ICA (spatial) with less options. This option can be accessed when you click on Setup ICA analysis button. Batch example is given in icatb/icatb_batch_files/batch_constrained_ica.m.
- Option is now provided to use an average mask in setup ICA analysis. Mask option can be accessed in “Setup-ICA defaults” menu.
- Some more dimensionality estimation options are provided in the Setup ICA analysis like:
    - MDL (FWHM): This option skips i.i.d sampling. You need to enter smoothness FWHM kernel used on the fMRI data.
    - Order estimated by entropy rate based methods (finite memory length and AR signal).
- Some more despike options are provided like despike based on smoothed timecourses as reference signal and median filtering. You can change these options in variable DESPIKE_OPTIONS in icatb_defaults.m.
- Batch option to do univariate tests directly is provided in the Mancovan toolbox. Options are provided to handle missing subjects at a particular voxel, frequency bin or FNC component pairs. Example templates are given in icatb/icatb_batch_files/input_mancovan_ttests.m.
- Option is now provided to use GIFTI data as input in SBM toolbox.
- Options are provided to merge separate ICA analyses in the Mancovan or dFNC along the subject dimension or component dimension (model order analysis given the same subjects). For more information, please see icatb/icatb_batch_files/input_dfnc.m file.
- DFNC related updates:
    - We added some features in the temporal dFNC toolbox.
        - Option is provided to do temporal variation FNC (Flor A. Espinoza et al., “Characterizing Whole Brain Temporal Variation of Functional Connectivity via Zero and First Order Derivatives of Sliding Window Correlations”, Front. Neurosci., 27 June 2019).
        - Temporal variation of functional network connectivity uses derivative of dFNC and searches for concurrent patterns in dFNC and its derivatives.
        - Average sliding window correlation is added (Victor M. Vergara et al., “An average sliding window correlation method for dynamic functional connectivity”, HBM, 2019).
        - Shared trajectory option is added (Ashkan Faghiri et al., “Weighted average of shared trajectory: A new estimator for dynamic functional connectivity efficiently estimates both rapid and slow changes over time”, J Neurosci Methods, 2020). Shared trajectory uses gradients to calculate weighted average of shared trajectory.
        - Model based dFNC is added (Ünal Sakoğlu et al., “A method for evaluating dynamic functional network connectivity and task-modulation: application to schizophrenia”, MAGMA, 2010). Task load function is computed at each window and correlated with the windowed correlations for each regressor.
    - Dynamic coherence toolbox (Maziar Yaesoubi et al., Dynamic coherence analysis of resting fMRI data to jointly capture state-based phase, frequency, and time-domain information, Neuroimage. 2015). Wavelet transform is used to compute dFNC in both time and frequency space.
    - Added an option to use ROI based dFNC. There are two options like ROI-ROI and ROI-voxel. ROI-ROI computes cross correlation between the averaged timecourses at each ROI. ROI-voxel computes cross-correlation between the average timecourse of given ROI and the rest of the brain at each voxel.
    - Windowless Functional connectivity (Maziar Yaesoubi et al., A window‐less approach for capturing time‐varying connectivity in fMRI data reveals the presence of states with variable rates of change, 2018). This approach calculates dFNC states as the outer product the bases of subspace estimated using K-SVD.
    - Spatial chronnectome toolbox is added. (Iraji, A. et al. (2019) 'The spatial chronnectome reveals a dynamic interplay between functional segregation and integration', Hum Brain Mapp, 40 (10), pp. 3058-3077). Spatial chronnectome captures voxel wise changes in the spatial patterns across time.
    - Spatial dynamics hierarchy toolbox (Iraji, A et al. (2019) 'Spatial dynamics within and between brain functional domains: A hierarchical approach to study time-varying brain function', Hum Brain Mapp, 40 (6), pp. 1969-1986). Spatial dynamics approach studies dynamic properties within the brain hierarchy.
 
***

## GroupICAT v4.0b (Feb 20, 2017):
- A subject outlier detection tool is now added as a "Generate Mask" utility in the GIFT toolbox. An average mask is generated and subjects below a certain correlation threshold are excluded from the analysis. At the end of the mask generation, a GIFT batch file is created which can be used to run the group ICA.
- An option is now provided to generate results summary in the Mancovan toolbox. This tool can be accessed using "display" button in the Mancovan toolbox. Univariate results are plotted in a separate figure for each significant covariate and a connectogram display (Rashid, B., et al. (2014), "Dynamic connectivity states estimated from resting fMRI Identify differences among Schizophrenia, bipolar disorder, and healthy control subjects", Frontiers in human neuroscience) is used to show the FNC plots.
- We now provide an option in the Mancovan toolbox to run univariate tests using the selected covariates bypassing the multivariate tests. This tool can be accessed using "Run analysis" button in the Mancovan toolbox.
- The "Group networks" tool is renamed to a more general "network summary" in the GIFT display tools. The network summary display uses the component network information and optional FNC information to generate composite orthogonal views (Damaraju, E., et al. (2014), "Dynamic functional connectivity analysis reveals transient states of dysconnectivity in schizophrenia", NeuroImage), composite rendered surfaces of brain, stacked orthogonal slices, FNC matrix viewer and connectogram FNC plot. Please see this HTMLreport for more information and example figures.
- We now provide an option to use the temporal design matrix information in the "Results Summary" button to compute R-square and one sample t-test on beta weights in the GIFT toolbox. This tool can also be accessed as temporal sorting under "Utilities" drop down box.
- The stand-alone image viewer tool is now enabled to select multiple component images which can be plotted independently or in a composite plot (montage, render or orthogonal slices).
- An option is now provided to export results to PDF or HTML file in the component viewer display tool.
- MOO-ICAR algorithm name is changed to GIG-ICA.

***

## GroupICAT v4.0a (May 02, 2015):
- Group ICA command line tool ("gica_cmd") is now added to the toolbox. Batch script can be run with minimal options from the MATLAB command line.
- Two more PCA methods are now integrated in the GIFT toolbox (early work based on: S. Rachakonda and V. D. Calhoun, "Efficient Data Reduction in Group ICA Of fMRI Data," in Proc. HBM, Seattle, WA, 2013, and in two papers currently under review)
    - Option to do PCA using Multi power iteration (MPOWIT) is now added. MPOWIT accelerates subspace iteration approach to extract dominant components from the data with very high accuracy in only a few iterations. MPOWIT can be run with data available in memory or by loading one data-set at a time.
    - Subsampled time PCA (STP) based on three data reduction method is now integrated in the toolbox. STP avoids whitening in the intermediate PCA step and PCA subspace is updated based on previous group estimates and new group entered. [note this addresses previous issues related to performance of three-step PCA]
- Parallel computing is now incorporated. Group ICA makes use of parallel computing toolbox to speed up the analysis stage. If parallel computing toolbox is not available, independent MATLAB sessions are used to speed up the process. The following tools are run in parallel:
    - Dimensionality Estimation
    - Subject level PCA
    - Stability analysis using ICASSO or Minimum spanning tree (MST)
    - Back reconstruction using spatial temporal regression or MOO-ICAR methods
    - Scaling components
    - Removing components
- Standalone image display tools like montage, orthogonal viewer, rendering and grouping components by network names are now integrated.
- Option is provided to summarize group ICA results using an HTML report.
- Spatial dynamic functional connectivity toolbox (sDFNC) is integrated. sDFNC toolbox is based on paper S. Ma, V. Calhoun, R. Phlypo and T. Adali, “Dynamic changes of spatial function network connectivity in healthy individuals and schizophrenia patients using independent vector analysis”, NeuroImage, 90 (2014), 196-206.
- Options are now provided to run both temporal and spatial ICA on pre-processed data directly.

***

## GroupICAT v3.0a (May 21, 2013):
- A new dynamic functional network connectivity (dFNC) toolbox is integrated within the GIFT toolbox. Please see E. Allen, E. Damaraju, S. M. Plis, E. Erhardt, T. Eichele, and V. D.Calhoun, "Tracking whole-brain connectivity dynamics in the resting state", Cereb Cortex, in press. The approach has also recently been validated using concurrent EEG/fMRI (see E. Allen, T. Eichele, L. Wu, and V. D. Calhoun, "EEG Signatures of Functional Connectivity States," in Human Brain Mapping, Seattle, WA, 2013).
- A new algorithm for independent vector analysis (IVA-GL) is integrated into the GIFT toolbox. IVA-GL uses multi-variate gaussian prior and laplacian prior to do source separation from the data. Please see M. Anderson, T. Adali, & X.-L. Li, "Joint Blind Source Separation of Multivariate Gaussian Sources: Algorithms and Performance Analysis," IEEE Trans. Signal Process., 2012, 60, 1672-1683) and T. Kim, H. T. Attias, S.-Y. Lee, & T.-W. Lee, "Blind Source Separation Exploiting Higher-Order Frequency Dependencies," IEEE Trans. Audio Speech Lang. Process., 2007, 15, 70-79.
- Multivariate Objective Optimization ICA with Reference (MOO-ICAR) is now integrated in the GIFT toolbox. MOO-ICAR uses a no data reduction approach and aggregate component maps from previous group ICA analysis as reference to estimate sources of interest for each subject. Please see Y. Du, Y. Fan, "Group information guided ICA for fMRI data analysis", NeuroImage 69: 157-197 (2013).
- Constrained ICA (Spatial) approach is updated to allow for the no data reduction approach similar to MOO-ICAR method above.
- PCA using eigen decomposition method is modified to compute the covariance matrix along the smallest dimension of the data.
- Statistics tool to compute T-test, ANOVA and Multiple Regression on subject ICA loading coefficients have now been added to the SBM toolbox.

***

## GroupICAT v2.0e (July 08, 2011):
- Mancovan toolbox is integrated in GIFT. Mancovan toolbox does multivariate tests on ICA timecourse spectral power, spatial map intensity and functional network connectivity to determine the significant covariates which will be used later in the univariate tests. Please see E. Allen, E. Erhardt, E. Damaraju, W. Gruner, J. Segall, R. Silva, M. Havlicek, S. Rachakonda, J. Fries, R. Kalyanam, A. Michael, J. Turner, T. Eichele, S. Adelsheim, A. Bryan, J. R. Bustillo, V. P. Clark, S. Feldstein, F. M. Filbey, C. Ford, K. Hutchison, R. Jung, K. A. Kiehl, P. Kodituwakku, Y. Komesu, A. R. Mayer, G. D. Pearlson, J. Phillips, J. Sadek, M. Stevens, U. Teuscher, R. J. Thoma, and V. D. Calhoun, "A baseline for the multivariate comparison of resting state networks," Frontiers in Human Neuroscience, vol. 1, p. 12, 2011
- SBM toolbox is added to do source based morphometry. Source based morphometry is a useful tool to study the gray matter differences between patients and controls. Please see below for references:
    - L. Xu, K. Groth, G. Pearlson, D. Schretlen, and V. Calhoun, "Source Based Morphometry: The Use of Independent Component Analysis to Identify Gray Matter Differences with Application to Schizophrenia," Hum Brain Mapp, vol. 30, pp. 711-724, 2009.
    - A. Caprihan, C. Abbott, J. Yamamoto, G. D. Pearlson, N. Bizzozero, J. Sui, and V. D. Calhoun, "Source-based morphometry analysis of group differences in fractional anisotropy in schizophrenia," Brain Connectivity, In Press.
- Options are provided in the Group ICA Toolbox to write only the necessary output components information which will be used later in the display.
- SPM8 volume functions are used to read and write image data.
- While using scaling timecourses option in GIFT, average of top 1% voxels is used instead of the maximum spatial intensity.
- Default mask used in the dimensionality estimation step is generated using all subjects in the analysis.

***

## GroupICAT v2.0d (March 31, 2010):
- Added data pre-processing options like intensity normalization, variance normalization and removing mean time-series. Intensity normalization is recommended for maximizing the reliability and replicability of the components. Please see abstract for more detalis.
- Added expectation maximization option to PCA computation (S.Roweis, "EM algorithms for PCA and sensible PCA", Advances in Neural Information Processing Systems. 1998). EM PCA has fewer memory constraints compared to covariance based PCA and is the preferred method for very large data-set analysis.
- All the subsequent analysis MAT files after PCA with the exception of ICASSO will be converted to single precision if you have selected single precision in the PCA options. This is recommended for maximizing the number of data sets given a fixed amount of RAM.
- GICA3 back reconstruction method was released as an update to GroupICAT v2.0c. GICA3 is the recommended method for reconstructing individual subject components with the most accurate spatial maps and timecourses. GICA3 has two desirable properties that the sum of the subject spatial maps is the aggregate spatial map and the product of the time courses and spatial maps estimate the data to the accuracy of the PCA's. We have performed extensive comparison with 3 different PCA/data-reduction approaches and 4 back-reconstruction approaches including spatio-temporal regression/dual regression. These approaches are included as options in GIFT. Please see abstract for more details.
- Component images sign in ICA is flipped based on the skewness measure of the distribution (previously it was flipped based upon the maximum voxel).
- Subject component image distributions are centered to zero by default when scaling component images. Centering is done based on the peak of the distribution.
- "icatb_mem_ica.m" script is updated to include all the data reduction strategies and will give a close estimate of how much RAM is required for all the analysis types.
- ICASSO can now be acessed from Setup ICA GUI.
- Dimensionality estimation source code ("icatb_estimate_dimension.m") developed by Leo Li is now available.
- Batch script doesn't read "numOfPC3" variable if only two data reduction steps are used.
- Bug was fixed in the back reconstruction code (GroupICAT v2.0c Updates, March 11, 2010) to handle "Constrained ICA (Spatial)" algorithm.
- Covariance options is replaced with PCA options and will be available when you select the PCA type in Setup ICA GUI.

***

## GroupICAT v2.0c (August 17, 2009):
- Enabled two data reduction steps in Setup ICA when the number of subjects is greater than 10.
- Two data reduction steps method is handled better in terms of memory usage for analyzing very large datasets.
- Added C-MEX files for computing eigen values of a symmetric matrix using packed storage scheme (this approach is slightly slower, but less memory intensive).
- Added spatial-temporal back reconstruction approach (this is an alternative approach to back-reconstruction and computes a spatial regression of the aggregate component images onto each timepoint of the single subject data and then computes a temporal regression of the single subject component timecourses onto each voxels timecourse). Overall results are quite similar to back-reconstruction using the PCA de-whitening matrix, however for the most accurate estimates of spatial maps we recommend GICA3.

***

## GroupICAT v2.0b (April 02, 2009):
- Single trial amplitudes utility based on Dr. Tom Eichele's work is now added to the GIFT toolbox.
- Added an option in the batch script to select the input data using regular expression pattern match. This can be used to get the directories that are highly nested.
- We remove the limitation to use MATLAB Statistics toolbox in order to compute statistics on the time courses.
- Added Multiple Regression design criteria in the Stats on time courses GUI.
- Percent variance calculation is added in the Utilities Section.

***

## GroupICAT v2.0a (April 11, 2008):
- We now provide EEGIFT toolbox for analyzing group ICA on EEG data (By Tom Eichele). EEGIFT contains options for importing data in .SET format from EEGLAB and visualization methods for viewing group ICA components. Both GIFT and EEGIFT are subsumed within GroupICAT v2.0a.
- Temporal sorting in GIFT is optimized. We load .MAT files for individual subject components instead of using time course images.
- Event average using deconvolution method (By Tom Eichele) is implemented.
- Event average utility in main figure window now has the options for selecting multiple regressors and components. Event average results will be written as .MAT files.
- Slider callback is optimized when very large number of time courses are plotted using "Split-timecourses" utility.

***

## GIFT v1.3d (Dec 18, 2007):
- GUI for doing statistical testing of time courses (beta weights) is included.
- Right-left text plotted during display is changed in this version to make it consistent with SPM convention (Neurological convention).
- Flip parameter for analyze images is stored in ICA parameter file and will give a warning message whenever flip parameter is changed during display.
- Statistics are done on component images and time courses even if the time points are different.
- Latest SPM updates for volume functions are installed.
- File selection window contains an option to enter a subset of Nifti files and an edit button to change the file selection.

***

## GIFT v1.3c (Jan 8, 2007):
- Constrained ICA (Spatial) algorithm developed by Qiu-Hua Lin is added to the GIFT toolbox.
- Default mask calculation is changed such that Boolean AND operation is performed on each data-set.
- PCA, ICA, Calibration MAT files contain only the voxels that are in the mask. Atleast 30% - 40% disk space will be saved.
- Both batch script and setup ICA GUI share the same code.
- Dimensionality estimation step is now batched.
- Display methods can now be accessed through a batch file.
- Error messages are reported with the line numbers on Matlab 7.

***

## GIFT v1.3b (April 21, 2006):
- Now writes 3D analyze images compatible with SPM2.
- Component images are detrended while converting to Z-scores.

***

## GIFT v1.3a (March 17, 2006) :
- Now reads functional data in Nifti or 4D Analyze format.
- New MDL Algorithm for estimating the number of components.
- Pre-compiled ICA algorithms: Erica, Simbec, Jade Opac, EVD and Amuse can now be run on Matlab 7.
- Option is provided to compress image files to zip format (to reduce the number of files and disk space).
- v1.3a reads the images from the previous version but older versions cannot read the images written using the new version.
- Correlation, regression, kurtosis and maximum voxel results are saved to a text file. The text file location will be printed to the command prompt.
- Manual is updated to include additional examples of temporal sorting and using output regression parameters, as well as statistical analysis of images using SPM2.
- Display GUI and setup ICA GUI are changed to minimize the selection process.
- Option is provide to calculate stats and event average under Utilities drop down box.

***

## GIFT v1.2d (December 5, 2005):
- Dimensionality estimation code is fixed to handle negative voxel dimensions and PCA is run on voxels that surpass the threshold.

***

## GIFTv1.2c with updates (November 7, 2005):
- Slices in mm are shown in Component Explorer and Composite viewer visualization methods.
- Option is provided in temporal sorting to automatically sort components or enter the regressor names through a text file.

***

## GIFTv1.2c with updates (October 3, 2005):
- Fixed bug for the component time courses that look flipped after calibration.
- Removing artifacts from the data is included.
- Temporal sorting with session specific regressors using one SPM2 design matrix is provided.
- Detrending of ICA time courses during scaling of components is included.

***

## GIFT v1.2c with the release date 5 July 2005:
- Data-sets can be analyzed with different number of images or scans but voxel dimensions should be the same.
- Error checking is done when SPM2 design matrix is loaded. Number of images of data-set is checked with the nscans field in SPM structure.
- Regressors specific to session can be selected during temporal sorting.
- Higher order detrending is provided when the components are sorted temporally.
- Multiple Regression step is optimized in sorting components.

***

## GIFT v1.2b with the release date 18 March 2005:
- Semi-blind ICA algorithm created by Vince Calhoun is added to the toolbox.
- Data reduction step is optimized in the ICA analysis step. Estimating components and the group ICA run faster than the previous version.
- New user interface is provided to select the reference functions while sorting the components.
- After the components are temporally sorted using Multiple Regression as sorting criteria, ICA Time courses can be adjusted by right clicking on the axis. ICA time courses are adjusted by removing the nuisance parameters and the regression coefficients other than the selected reference function.
- In the event related average option is provided to select the reference function.
- In batch scripts option is provided to specify one design matrix for all subjects. Components can be visualized using the Display GUI.
- Fixed bug in batch scripts when subject folder names with variable length are specified.

***

## GIFT v1.2a with the release date 26 November 2004:
- New user interface for Setup ICA is provided. ICA parameters and algorithms can be easily switched. Help button is included adjacent to each parameter.
- Event related average for the ICA time course is calculated based on the onsets of the given SPM model or design matrix.
- Time course window for multiple subjects and sessions is shown in a new figure window with a scroll bar.
- New color maps for the composite viewer are provided. A maximum of five different color bars are plotted.
- Batch script with two sample text files is provided. The explanation of the parameters is given in the help manual.
- Interactive file selection window is updated to show the number of selected files and directory history.
- Regression parameters or the coefficients are written to .mat and .txt files when the components are sorted with Multiple Regression sorting criteria.

***

## GIFT v1.1d with the release date 14 October 2004 (Released for the ICA Class at the Olin Neuropsychiatry Research Center):
- Estimating number of independent components from the fMRI data is included. The number of components estimated is shown to the user before selecting the number of principal components.
- HTML Help Manual is provided. When the help button is pressed HTML Help is opened in the default browser.
- Maximum Voxel sorting criteria is added to sort components spatially based on the spatial template selected.
- Optimal ICA algorithm created by Baoming Hong and Vince Calhoun is added.
- Interactive file selection window is added instead of using the spm_get function.

***

## GIFT v1.01d with the release date 13 September 2004:
- Components can be spatially sorted by using a template image which contains regions of interest. Sample templates are provided in the folder 'icatb_templates' with names 'LeftTemplate.img' and 'RightTemplate.img'.
- All the analysis information is stored in a log file which ends with '_results.log'. This file gets appended every time when you run the analysis with the same prefix for the output files.
- The parameters involved in the ICA and PCA step are shown to the user when the subject files are already selected.
- Scroll bar is provided for the edit and pop up controls in the input dialog box for selecting the ICA options.
- Any error occurring during setting up the analysis, running the analysis or displaying the results is shown to the user.
- Fixed bug with correlation sorting criteria - When 'temporal' is selected in 'Select Sorting Type' option and 'Select All Subjects and Sessions' is selected in 'What do you want to sort?' option.
- Horizontal scroll bar in all the dialog boxes is turned off.

***

## GIFT v1.01c with the release date 20 August 2004:
- Number of partitions option in Setup ICA Analysis is turned off.
- Normalize model time course checkbox option removed for the components that are not sorted in the figure window which will be displayed when clicked on the time course window.
- Defaults in the 'icatb_defaults.m' file are applied to display GUI window. Option is provided to the user to change the defaults. Time courses can be smoothed by replacing the parameter 'SMOOTHPARA' from 'no' to 'yes' and the value can be changed by giving different values to 'SMOOTHINGVALUE' parameter in 'icatb_defaults.m' file.Four options for detrending time courses are provided in case of sorting with the Multiple Regression sorting criteria.
- The 'DETRENDNUMBER' parameter in 'icatb_defaults.m' file can be changed from '0' to '3' depending on the type of detrending you want to do. Comments are included in the 'icatb_defaults.m' file which explains what 'DETRENDNUMBER' means.

***

## GIFT  v1.01b with the release date 28 July 2004:
- Fixed bug for the components to be sorted when you use combination of "No design matrix" in the Setup ICA Analysis and "select model for every subject" in Sort Components GUI. New dialog boxes for showing the directions about the toolbox.
- New input dialog box for ICA algorithms which can incorporate any number of inputs from the user.
- Detrend is done prior to concatenation of the time courses in sorting components.
- Line fit is shown along with the model and ICA time courses for Multiple Regression sorting criteria. Help button which explains how to use Group ICA Toolbox.
- Status bar in the run analysis button which shows how much percentage of the analysis is done.ICA algorithms like Simbec, Evd, Jade Opac and Amuse are included in compiled version.
- Option for selecting one regressor or multiple regressors is removed in sorting components GUI.
