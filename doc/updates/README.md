# **GIFT Updates**

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
