GroupICATv4.0b Updates (Oct 19, 2018):

1. Error message "Undefined function 'mtimes' for input arguments of type 'cell" is fixed when running SBM toolbox.

2. Option is provided to use temporary directory when unzipping gzip files in gzip nifti reader. 

3. GIFT Nipype plugin is updated to use in the latest version of Nipype.

4. GIG-ICA function is optimized.

The following files are updated:

1. icatb/icatb_defaults.m
2. icatb/icatb_analysis_functions/icatb_dataReduction.m
3. icatb/icatb_analysis_functions/icatb_algorithms/icatb_gigicar.m
4. icatb/icatb_io_data_functions/icatb_read_gzip_nii.m
5. icatb/nipype-0.10.0/nipype/interfaces/gift/base.py


GroupICATv4.0b Updates (Feb 09, 2018):

PCA can now be done without loading the entire data when using the SBM toolbox. Please set PCA options appropriately for large data. Function icatb_dataReduction is changed.

GroupICATv4.0b Updates (Jan 23, 2018):

Sign of the components is changed based on the skewness when using SBM toolbox. File icatb/icatb_analysis_functions/icatb_calculateICA is updated.


GroupICATv4.0b Updates (Dec 04, 2017):

Index exceeds dimensions error is fixed when computing FNC correlations in the mancovan toolbox. Function icatb/icatb_helper_functions/icatb_compute_fnc_corr.m is updated.

GroupICATv4.0b Updates (Dec 01, 2017):

An option is now provided to do MDL estimation in mancovan toolbox:

The following files are updated:

1. icatb/icatb_mancovan_files/icatb_run_mancovan.m
2. icatb/icatb_mancovan_files/icatb_setup_mancovan.m
3. icatb/icatb_mancovan_files/icatb_mancovan_batch.m

GroupICATv4.0b Updates (Nov 14, 2017):

1. GIFT interface for Nipype is now added (http://nipype.readthedocs.io/en/latest/). For more information on how to use GIFT on Nipype, please follow readme.txt file in icatb/nipype-0.10.0/nipype/interfaces/gift
2. Option is now provided to use lag information when computing FNC features in the Mancovan toolbox.
3. Nifti Gzip memory option is set to handle very large files with less heap space usage.
4. Option is provided to generate HTML/PDF reports using a batch file.
5. Network summary is updated to use parameter file in addition to nifti/FNC correlation files.
6. Create mask is updated to use Nifti Gzip when running ICA in parallel mode.
7. Connectogram function is updated to use outer radius of circle to plot the spatial maps closer to the circle. 
8. By default, SVD is run when using the first PCA step to both maximize performance and reduce memory usage.
9. Option is provided to do connectogram plots on the centroids from standard dFNC.
10. Options are provided to bypass i.i.d sampling and use smoothness kernel of the data when doing dimensionality estimation. Please see variable DIM_ESTIMATION_OPTS in icatb_defaults.m
11. Option is now provided to run mancova on the network domain averaged static FNC correlations. Please see variable MANCOVA_DEFAULTS.fnc.domain_average to run mancova on domain averaged FNC correlations.

The following files are updated/added:

1. icatb/icatb_defaults.m
2. icatb/post_process_dfnc.m
3. icatb/gica_cmd.m
4. icatb/icatb_analysis_functions/icatb_runAnalysis.m
5. icatb/icatb_batch_files/icatb_batch_file_run.m
6. icatb/icatb_batch_files/icatb_read_batch_file.m
7. icatb/icatb_batch_files/Input_spatial_ica.m
8. icatb/icatb_display_functions/icatb_network_summary_gui.m
9. icatb/icatb_display_functions/icatb_plot_connectogram.m
10. icatb/icatb_display_functions/icatb_resizeImage.m
11. icatb/icatb_display_functions/icatb_disp_temp_regress_results.m
12. icatb/icatb_display_functions/icatb_component_viewer.m
13. icatb/icatb_helper_functions/icatb_generateMask.m
14. icatb/icatb_helper_functions/icatb_dfnc_options.m
15. icatb/icatb_helper_functions/icatb_gica_html_report.m
16. icatb/icatb_helper_functions/icatb_optimal_clusters.m
17. icatb/icatb_helper_functions/icatb_loadAndInterpTC.m
18. icatb/icatb_helper_functions/icatb_run_dfnc.m
19. icatb/icatb_helper_functions/icatb_dfnc_results.m
20. icatb/icatb_helper_functions/icatb_model_order_prediction.m
21. icatb/icatb_helper_functions/icatb_postprocess_timecourses.m
22. icatb/icatb_helper_functions/icatb_sbm_html_report.m
23. icatb/icatb_helper_functions/icatb_interp_data.m
24. icatb/icatb_helper_functions/icatb_display_dfnc.m
25. icatb/icatb_helper_functions/icatb_report_generator.m
26. icatb/icatb_helper_functions/icatb_compute_fnc_corr.m
27. icatb/icatb_helper_functions/icatb_filt_data.m
28. icatb/icatb_io_data_functions/icatb_OptionsWindow.m
29. icatb/icatb_io_data_functions/icatb_read_gzip_nii.m
30. icatb/icatb_mancovan_files/icatb_display_mancovan.m
31. icatb/icatb_mancovan_files/icatb_fitggmix.m
32. icatb/icatb_mancovan_files/icatb_mancovan_feature_options.m
33. icatb/icatb_mancovan_files/icatb_mysplinefun.m
34. icatb/icatb_mancovan_files/icatb_plot_mancova_features.m
35. icatb/icatb_mancovan_files/icatb_plot_mult_mancovan.m
36. icatb/icatb_mancovan_files/icatb_plot_univariate_results.m
37. icatb/icatb_mancovan_files/icatb_run_mancovan.m
38. icatb/icatb_mancovan_files/icatb_setup_mancovan.m
39. icatb/icatb_parallel_files/icatb_parCreateMask.m
40. icatb/icatb_parallel_files/icatb_par_postprocess_timecourses.m
41. icatb/icatb_helper_functions/icatb_estimateCompCallback.m
42. icatb/icatb_estimate_dimenison.m
43. icatb/nipype-0.10.0 directory is added.
44. icatb/mancovan_toolbox.m and icatb/mancovan_toolbox.fig.
45. icatb/icatb_helper_functions/icatb_domain_avg_fnc.m 
46. icatb/icatb_helper_functions/icatb_low2Mat.m
47. icatb/icatb_mancovan_files/icatb_mancovan_batch.m


GroupICATv4.0b Updates (May 01, 2017):

spm mexmaci64 files are updated to work on R2017a. 

GroupICATv4.0b Updates (April 28, 2017):

PCA paper reference S. Rachakonda, R. F. Silva, J. Liu and V. D. Calhoun, "Memory Efficient PCA Methods for Large Group ICA", Frontiers in Neuroscience, 2016 is mentioned in the 
icatb/icatb_analysis_functions/icatb_calculate_pca.m and icatb/icatb_parallel_files/icatb_parCalculatePCA.m files.

GroupICATv4.0b Updates (April 17, 2017):

icatb/icatb_display_functions/icatb_plot_connectogram.m file is updated to include:

	1. By default, only 12 network colors are provided instead of 14. Line number 244 is fixed.
	2. An option is provided to set connectogram figure size in icatb_defaults.m file. Please see CONNECTOGRAM_FIG_POS variable.
	3. Optional flag is provided not to reorder the connectivity matrix when all the components are selected.

GroupICATv4.0b Updates (April 13, 2017):

1. icatb/icatb_io_data_functions/icatb_read_gzip_nii.m is fixed to load data correctly when slices are entered as an input in the function.
2. Components (and its labels) were not updated in the connectogram after removing zeros in the FNC matrix (April 11th update). icatb/icatb_display_functions/icatb_plot_connectogram function 
is now fixed. 


GroupICATv4.0b Updates (April 12, 2017):

Network color is set to a default value if uisetcolor window is closed without selecting a value. Function icatb/icatb_display_functions/icatb_plot_connectogram.m is updated.

GroupICATv4.0b Updates (April 11, 2017):

1. GIFT now reads NIFTI Gzip (*.nii.gz) format. For very large files (multiband data) it is recommended to set NIFTI_GZ = 1 in icatb_defaults.m (files will be un-archived) for reading 
Gzip files faster.
2. Despike utility computational performance (speed) is improved.
3. We now provide various options to estimate clusters using elbow, BIC, AIC and Dunns index. Optimal clusters is also provided as a stand alone tool.
4. Connectogram is now added as a stand alone tool in GIFT -> Tools -> Display Tools -> Misc menu. Options to use rendered images, user defined RGB images and component labels are provided 
in the connectogram function. 
5. Generate mask utility is now added in the SBM toolbox. To access the tool, use SBM -> Utilities -> Generate Mask.
6. Deprecated legend parameter is removed in icatb_orthoViewer.m.

The following files are added or updated:	

	1. icatb/dfnc_toolbox.fig      
	2. icatb/dfnc_toolbox.m        
	3. icatb/gift.fig              
	4. icatb/icatb_defaults.m      
	5. icatb/icatb_displayGUI.m    
	6. icatb/icatb_setup_analysis.m
	7. icatb/sbm.fig   
	8. icatb/icatb_analysis_functions/icatb_parameterInitialization.m
	9. icatb/icatb_analysis_functions/icatb_runAnalysis.m   
	10. icatb/icatb_batch_files/icatb_read_batch_file.m
	11. icatb/icatb_display_functions/icatb_display_composite.m
	12. icatb/icatb_display_functions/icatb_orthoViewer.m
	13. icatb/icatb_display_functions/icatb_plot_connectogram.m
	14. icatb/icatb_helper_functions/icatb_createMask.m         
	15. icatb/icatb_helper_functions/icatb_despike.m            
	16. icatb/icatb_helper_functions/icatb_despike_tc.m         
	17. icatb/icatb_helper_functions/icatb_generateMask.m       
	18. icatb/icatb_helper_functions/icatb_get_countTimePoints.m
	19. icatb/icatb_helper_functions/icatb_loadAndInterpTC.m    
	20. icatb/icatb_helper_functions/icatb_optimal_clusters.m   
	21. icatb/icatb_helper_functions/icatb_post_process_dfnc.m  
	22. icatb/icatb_helper_functions/icatb_read_hdr.m           
	23. icatb/icatb_helper_functions/icatb_rename_4d_file.m     
	24. icatb/icatb_helper_functions/icatb_returnHInfo.m        
	25. icatb/icatb_helper_functions/icatb_update_mask.m        
	26. icatb/icatb_helper_functions/icatb_utilities.m  
	27. icatb/icatb_io_data_functions/icatb_dataSelection.m
	28. icatb/icatb_io_data_functions/icatb_gz.jar         
	29. icatb/icatb_io_data_functions/icatb_inputdlg2.m    
	30. icatb/icatb_io_data_functions/icatb_read_data.m    
	31. icatb/icatb_io_data_functions/icatb_read_gzip_nii.m
	32. icatb/icatb_mancovan_files/icatb_mysplinefun.m 
	33. icatb/icatb_mancovan_files/icatb_run_mancovan.m
	

GroupICATv4.0b Updates (March 01, 2017):

icatb/icatb_mancovan_files/icatb_run_mancovan.m file is fixed to handle NaNs in the nuisance covariates.