GroupICATv4.0c Updates (Nov 18, 2020):

1. BIDS data support is now added in the GIFT toolbox. 
2. Added an option to use second level IVA analysis (IVA-G) in the algorithm IVA-GL.
3. Option is now added to use tall array when computing K-means in spatial chronnectome toolbox.

The following files are updated or added:

        icatb\icatb_setup_analysis.m                                         
        icatb\post_process_spatial_chronnectome.fig                          
        icatb\post_process_spatial_chronnectome.m                            
        icatb\setup_reference_ica.m                                          
        icatb\icatb_analysis_functions\icatb_dataReduction.m                 
        icatb\icatb_analysis_functions\icatb_icaAlgorithm.m                  
        icatb\icatb_analysis_functions\icatb_icaOptions.m                    
        icatb\icatb_analysis_functions\icatb_parameterInitialization.m       
        icatb\icatb_analysis_functions\icatb_runAnalysis.m                   
        icatb\icatb_batch_files\Input_spatial_ica.m                          
        icatb\icatb_batch_files\Input_spatial_ica_bids.m                     
        icatb\icatb_batch_files\icatb_batch_file_run.m                       
        icatb\icatb_batch_files\icatb_read_batch_file.m                      
        icatb\icatb_helper_functions\icatb_compute_dfnc.m                    
        icatb\icatb_helper_functions\icatb_copy_ica_results_to_bids.m        
        icatb\icatb_helper_functions\icatb_generateMask.m                    
        icatb\icatb_helper_functions\icatb_kmeans_clustering.m               
        icatb\icatb_helper_functions\icatb_merge_analyses.m                  
        icatb\icatb_helper_functions\icatb_parseBIDS.m                       
        icatb\icatb_helper_functions\icatb_postprocess_sdh.m                 
        icatb\icatb_helper_functions\icatb_postprocess_spatial_chronnectome.m
        icatb\icatb_helper_functions\icatb_postprocess_timecourses.m         
        icatb\icatb_helper_functions\icatb_run_spatial_chronnectome.m        
        icatb\icatb_helper_functions\icatb_select_bids_params.m              
        icatb\icatb_helper_functions\icatb_update_mask.m                     
        icatb\icatb_io_data_functions\icatb_dataSelection.m                  
        icatb\icatb_parallel_files\icatb_par_postprocess_timecourses.m       
        icatb\icatb_spm_files\icatb_spm_BIDS.m                               
        icatb\icatb_spm_files\icatb_spm_jsonread.mexa64                      
        icatb\icatb_spm_files\icatb_spm_jsonread.mexmaci64                   
        icatb\icatb_spm_files\icatb_spm_jsonread.mexw32                      
        icatb\icatb_spm_files\icatb_spm_jsonread.mexw64                      
        icatb\icatb_spm_files\icatb_spm_load.m                               
        icatb\icatb_spm_files\icatb_write_nifti_data.m                       
        icatb\nipype-0.10.0\nipype\interfaces\gift\model.py                  
        icatb\nipype-0.10.0\nipype\interfaces\gift\readme.txt
	icatb\toolbox\noisecloud\noisecloud_spm_coregister.m                


