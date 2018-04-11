function inputText = icatb_define_parameters(varargin)
% The input parameters are defined below. read_from_file field states that
% the field can be read from input file if provided.
%
% Inputs:
% inputFile - Input file required for batch analysis
%
% Output:
% inputText - Structure containing necessary parameter information


% Load defaults
icatb_defaults;
global PREPROC_DEFAULT;
global PCA_DEFAULT;
global BACKRECON_DEFAULT;
global SCALE_DEFAULT;
global GROUP_PCA_DEFAULT;


groupPCAOpts = char('Subject Specific', 'Grand Mean');

% Return modality type
modalityType = icatb_get_modality;

group_pca_default = matchString(GROUP_PCA_DEFAULT, groupPCAOpts, 1);
preproc_default = matchString(PREPROC_DEFAULT, icatb_preproc_data, 1);
pca_default = matchString(PCA_DEFAULT, icatb_pca_options, 1);

%% Remove GICA3 and STR from the list
bd = BACKRECON_DEFAULT;
if (ischar(bd) && strcmpi(bd, 'str'))
    bd = 'spatial-temporal regression';
end
backReconOptions = cellstr(icatb_backReconOptions);
if (isnumeric(bd))
    if (bd > length(backReconOptions))
        error(['Value for BACKRECON_DEFAULT exceeded the no. of backreconstruction (', num2str(length(backReconOptions)), ') options available']);
    end
    bd = deblank(backReconOptions{bd});
end
bd = lower(bd);
bad_inds = strcmpi(backReconOptions, 'regular') | strcmpi(backReconOptions, 'gica3');
backReconOptions(bad_inds) = [];
backrecon_default = matchString(bd, backReconOptions, 2);

sd = SCALE_DEFAULT;

if (isempty(sd))
    sd = 0;
end

if (ischar(sd))
    if (strcmpi(sd, 'no'))
        sd = 'no scaling';
    elseif (findstr(sd, 'norm ics'))
        if (strcmpi(modalityType, 'fmri'))
            sd = 'scaling in timecourses';
        else
            sd = 'scaling in topographies';
        end
    elseif (findstr(sd, 'joint scaling'))
        if (strcmpi(modalityType, 'fmri'))
            sd = 'scaling in maps and timecourses';
        else
            sd = 'scaling in timecourses and topographies';
        end
    elseif (strcmpi(sd, 'calibrate'))
        if (strcmpi(modalityType, 'fmri'))
            sd = 'scale to original data(%)';
        else
            sd = 'scale to original data';
        end
    end
else
    sd = sd + 1;
end

scaleOptions = icatb_scaleICA;

if (strcmpi(modalityType, 'fmri'))
    scale_default = matchString(sd, scaleOptions, 4);
else
    scale_default = matchString(sd, scaleOptions, 3);
    newScaleOptions = char('No Scaling', 'Z-scores');
    scale_default = strmatch(lower(deblank(scaleOptions(scale_default, :))), lower(newScaleOptions), 'exact');
    if (isempty(scale_default))
        scale_default = 1;
    end
    scaleOptions = newScaleOptions;
end

%%% Define the parameters for the set up ICA analysis %%%%%
% Parameter 1
numParameters = 1;

inputText(numParameters).promptString = 'Enter Name(Prefix) Of Output Files';
inputText(numParameters).answerString = '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'prefix';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Output Files', 'string', 'All the output files will be preprended with this prefix.');

numParameters = numParameters + 1;

inputText(numParameters).promptString = ['Have You Selected The ', modalityType, ' Data Files?'];
inputText(numParameters).answerString = 'Select';
inputText(numParameters).uiType = 'pushbutton';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'files';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
if (~strcmpi(modalityType, 'smri'))
    inputText(numParameters).help = struct('title', 'Data', 'string', ['There are two ways to select the data. Option one requires a separate directory for each subject and all subject directories are', ...
        ' stored in one root directory and subjects have the same file pattern. Otherwise use option 2. After the data is selected, *Subject.mat file is created and the component numbers are by default set to 20.']);
else
    inputText(numParameters).help = struct('title', 'Data', 'string', 'Select all subject images. After the data is selected, *Subject.mat file is created and the component numbers are by default set to 20.');
end


if (~strcmpi(modalityType, 'smri'))
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString =  'Select Type Of Data Pre-processing';
    inputText(numParameters).answerString = char(icatb_preproc_data);
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'preproc_type';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = preproc_default;
    inputText(numParameters).flag = 'scalar';
    inputText(numParameters).help = struct('title', 'Pre-processing', 'string', char('Data is pre-processed prior to the first data reduction. Options are discussed below:', '1) Remove Mean Per Timepoint - At each time point, image mean is removed.', ...
        '2) Remove Mean Per Voxel - Time-series mean is removed at each voxel', '3) Intensity Normalization - At each voxel, time-series is scaled to have a mean of 100. Since the data is already scaled to percent signal change, there is no need to scale the components.', ...
        '4) Variance Normalization - At each voxel, time-series is linearly detrended and converted to Z-scores.'));
end

if (~strcmpi(modalityType, 'eeg'))
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString =  'What Mask Do You Want To Use?';
    inputText(numParameters).answerString =  char('Default Mask', 'Select Mask');
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'maskFile';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    inputText(numParameters).flag = 'scalar';
    
    if (strcmpi(modalityType, 'fmri'))
        inputText(numParameters).help = struct('title', 'Mask', 'string', char('1) Default Mask - Mask is calculated using all the files for subjects and sessions or only the first file for each subject and session depending upon the variable DEFAULT_MASK_OPTION value in defaults. Boolean AND operation is done to include the voxels that surpass the mean of each subject''s session.', ...
            '2) Select Mask - Masks must be in Analyze or Nifti format.'));
    else
        inputText(numParameters).help = struct('title', 'Mask', 'string', char('1) Default Mask - Default mask includes voxels greater than or equal to 1% of mean.', '2) Select Mask - Masks must be in Analyze or Nifti format.'));
    end
    
    
end


numParameters = numParameters + 1;

inputText(numParameters).promptString =  'Select Type Of PCA';
inputText(numParameters).answerString = char(icatb_pca_options);
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'pcaType';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = pca_default;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).userdata = [];
inputText(numParameters).help = struct('title', 'PCA', 'string', char('There are five options like Standard, Expectation Maximization, SVD, MPOWIT and STP.', ...
    '1) Standard - Eigen value decomposition method is used to determine dominant components of interest', ...
    '2) Expectation Maximization - Expectation maximization method involves expectation and maximization steps to determine PCA components. It has fewer memory requirements but can be very slow if EM PCA is solved by loading data-set at a time. It is preferred to use MPOWIT method for faster convergence.', ...
    '3) SVD - Singular value decomposition', ...
    '4) MPOWIT - Multi power iteration method accelerates subspace iteration method to determine dominant components. MPOWIT converges in only a few iterations and is very useful when large data needs to be analyzed', ...
    '5) STP - Subsampled time pca method is a variation of 3 step pca method which involves grouping subjects into groups and retaining intermediate components which are usually higher than the final number of components estimated'));

if (~strcmpi(modalityType, 'smri'))
    
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString =  'Select Type Of Group PCA';
    inputText(numParameters).answerString = groupPCAOpts;
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'group_pca_type';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = group_pca_default;
    inputText(numParameters).flag = 'scalar';
    if strcmpi(modalityType, 'fmri')
        inputText(numParameters).help = struct('title', 'Group PCA Type', 'string', char('1.) Subject Specific - PCA is done on each dataset before doing group PCA.', ...
            '2.) Grand Mean - Each dataset is projected on to the eigen space of the mean of all datasets before doing group PCA. Please make sure that you have selected equal no. of timepoints between datasets.'));
    else
        inputText(numParameters).help = struct('title', 'Group PCA Type', 'string', char('1.) Subject Specific - PCA is done on each dataset before doing group PCA.', ...
            '2.) Grand Mean - Each dataset is projected on to the eigen space of the mean of all datasets before doing group PCA. Please make sure that you have selected equal no. of electrodes between datasets.'));
    end
    
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString =  'Select The Backreconstruction Type';
    inputText(numParameters).answerString =  backReconOptions;
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'backReconType';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = backrecon_default;
    inputText(numParameters).flag = 'scalar';
    
    if strcmpi(modalityType, 'fmri')
        inputText(numParameters).help = struct('title', 'Backreconstruction', 'string', char('Options available are GICA, Regular (GICA2), GICA3, Spatial-temporal regression and GIG-ICA. GICA2 and GICA3 are not not shown in the GUI but can be called in the batch script. GICA is a more robust tool to back reconstruct components when compared to GICA2 and GICA3 for low model order.', '', ...
            '1) GICA - GICA uses PCA whitening and dewhitening information to back reconstruct the component maps and timecourses.', ...
            '2) Spatial-temporal Regression - Back reconstruction is done using a two step multiple regression. In the first step, aggregate component spatial maps are used as basis functions and projected on to the subject''s data resulting in subject component time courses. In the second step, subject component time courses are used as basis functions and projected on to the subject''s data resulting in component spatial maps for that subject', ...
            '3) GIG-ICA - Aggregate components from group ica step are used to reconstruct individual subject components using GIG-ICA algorithm.', ...
            '', 'Note:', '1) GICA, GICA2 timecourses are similar to the timecourses obtained using Spatial-temporal Regression.', '2) Spatial maps obtained using GICA2 are exactly equal to the GICA3 method.', '3) All the back reconstruction methods give the same spatial maps and timecourses for one single subject single session analysis.', '4) GICA and Spatial-temporal Regression component timecourses are equivalent when 100% variance is retained in the first step PCA.' ...
            ));
    else
        inputText(numParameters).help = struct('title', 'Backreconstruction', 'string', char('Options available are GICA, Regular (GICA2), GICA3, Spatial-temporal regression and GIG-ICA. GICA2 and GICA3 are not not shown in the GUI but can be called in the batch script. GICA is a more robust tool to back reconstruct components when compared to GICA2 and GICA3 for low model order.', '', ...
            '1) GICA - GICA uses PCA whitening and dewhitening information to back reconstruct the component topographies and timecourses.', ...
            '2) Spatial-temporal Regression - Back reconstruction is done using a two step multiple regression. In the first step, aggregate component timecourses are used as basis functions and projected on to the subject''s data resulting in subject component topographies. In the second step, subject component topographies are used as basis functions and projected on to the subject''s data resulting in component timecourses for that subject', ...
            '3) GIG-ICA - Aggregate components from group ica step are used to reconstruct individual subject components using GIG-ICA algorithm.', ...
            '', 'Note:', '1) GICA, GICA2 topographies are similar to the timecourses obtained using Spatial-temporal Regression.', '2) Timecourses obtained using GICA2 are exactly equal to the GICA3 method.', '3) All the back reconstruction methods give the same topographies and timecourses for one single subject single session analysis.', '4) GICA and Spatial-temporal Regression component topographies are equivalent when 100% variance is retained in the first step PCA.' ...
            ));
    end
    
end

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Do You Want To Scale The Results?';
inputText(numParameters).answerString = scaleOptions;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'scaleType';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = scale_default;
inputText(numParameters).flag = 'scalar';

if strcmpi(modalityType, 'fmri')
    
    inputText(numParameters).help = struct('title', 'Scaling', 'string', char('1) Scale to Original Data(%) - Components are scaled to original data units to represent percent signal change.', ...
        '2) Z-scores - Components are scaled to z-scores.', ...
        '3) Scaling in Timecourses - Normalize spatial maps using maximum value (not absolute value) and apply it to timecourses.', ...
        '4) Scaling in Maps and Timecourses - Spatial maps are scaled using the standard deviation of timecourses and timecourses are scaled using the maximum spatial intensity value.'));
    
elseif (strcmpi(modalityType, 'smri'))
    
    inputText(numParameters).help = struct('title', 'Scaling', 'string', 'Z-scores - Components are scaled to z-scores.');
    
else
    
    inputText(numParameters).help = struct('title', 'Scaling', 'string', char('1) Scale to Original Data - Components are scaled to original data units.', ...
        '2) Z-scores - Components are scaled to z-scores.', ...
        '3) Scaling in Topographies - Normalize timecourses using maximum amplitude of timecourses (not absolute value) and apply it to topographies.', ...
        '4) Scaling in Timecourses and Topographies - Timecourses are scaled using the standard deviation of topographies and topographies are scaled using the maximum amplitude of timecourses.'));
    
end



numParameters = numParameters + 1;

% Open the input dialog if the Infomax, Fast ICA, Optimal ICA
inputText(numParameters).promptString = 'Which Algorithm Do You Want To Use?';
inputText(numParameters).answerString = icatb_icaAlgorithm;
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'algoType';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';

inputText(numParameters).help = struct('title', 'ICA Algorithm', 'string', 'ICA detects the independent sources present in the data');


if ~strcmpi(modalityType, 'eeg')
    
    numParameters = numParameters + 1;
    
    % Define the callback for this popup
    inputText(numParameters).promptString =  'Do You Want To Estimate The Number Of Independent Components?';
    inputText(numParameters).answerString =  char('No', 'Yes');
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'numeric';
    inputText(numParameters).tag = 'estimate_components';
    inputText(numParameters).enable = 'off';
    inputText(numParameters).value = 1;
    inputText(numParameters).flag = 'scalar';
    inputText(numParameters).help = struct('title', 'Dimensionality Estimation', 'string', ['Components are estimated from the fMRI data using the MDL criteria.', ...
        'Mean of the estimated components for all data-sets is taken when you use all subjects to estimate the components and the statistics are performed. After the statistics are calculated, all the information with a Plot button will be shown in a dialog box. The mean of MDL plot for all data-sets will be shown when you click the Plot button.', ...
        'When a particular subject''s session is selected the independent components are estimated and the plot of MDL is shown when clicked on the Plot button in the dialog box.']);
    
    
end


isHandleVisible = 'on';
if strcmpi(modalityType, 'smri')
    isHandleVisible = 'off';
end


if (strcmpi(modalityType, 'fmri'))
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString =  'Select Group ICA Type';
    inputText(numParameters).answerString =  char('Spatial', 'Temporal');
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'group_ica_type';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    inputText(numParameters).flag = 'scalar';
    inputText(numParameters).help = struct('title', 'Group ICA Type', 'string', char('1) Spatial ICA - Independent components are determined by maximizing independence in space.', ...
        '2) Temporal ICA - Independent components are determined by maximizing independence in time.'));
end


numParameters = numParameters + 1;

inputText(numParameters).promptString =  'How Many Data Reduction(PCA) Steps Do You Want To Run?';
inputText(numParameters).answerString =  char('2', '1');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'numReductionSteps';
inputText(numParameters).enable = 'off';
inputText(numParameters).visible = isHandleVisible;
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Data Reduction Steps', 'string', 'The number of times you want to do PCA. The number of data reduction steps depends on the number of data sets. For single subject single session only one PCA is automatically done. Two data reductions is preferred for multiple subjects and sessions.');

numParameters = numParameters + 1;

inputText(numParameters).promptString =  'Number Of PC (Step 1)';
inputText(numParameters).answerString =  '0';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'numOfPC1';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).visible = isHandleVisible;
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'PC1', 'string', 'No. of components desired in the first PCA step.');

numParameters = numParameters + 1;

inputText(numParameters).promptString =  'Number Of PC (Step 2)';
inputText(numParameters).answerString =  '0';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'numOfPC2';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).visible = isHandleVisible;
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'PC2', 'string', 'No. of components desired in the second PCA step.');

numParameters = numParameters + 1;

inputText(numParameters).promptString =  'Number Of IC';
inputText(numParameters).answerString =  '';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).dataType = 'numeric';
inputText(numParameters).tag = 'numComp';
inputText(numParameters).enable = 'inactive';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'IC', 'string', 'No. of independent components. This is the same as the no. of components selected in the last PCA step.');

if (~strcmpi(modalityType, 'smri'))
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString =  'Do you want to autofill data reduction values?';
    inputText(numParameters).answerString =  char('No', 'Yes');
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'autofill';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    inputText(numParameters).flag = 'scalar';
    inputText(numParameters).help = struct('title', 'Autofill', 'string', 'Numbers for principal components are auto-filled to 20 when the data is selected.');
end


numParameters = numParameters + 1;

inputText(numParameters).promptString =  'Select stability analysis type';
inputText(numParameters).answerString =  char('Regular', 'ICASSO', 'MST');
inputText(numParameters).uiType = 'popup';
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'which_analysis';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;
inputText(numParameters).flag = 'scalar';
inputText(numParameters).help = struct('title', 'Stability Analysis Type', 'string', char('1) Regular - ICA is run only once on the data.', ...
    '2) ICASSO toolbox is used to determine the reliability of ICA algorithm. ICA algorithm is run several times to determine the algorithmic reliability or stability. Reliable estimates correspond to tight clusters and un-reliable ones do not point to any cluster.', ...
    '3) MST - Uses minimum spanning tree to determine stable ICA run estimates.'));


if (~strcmpi(modalityType, 'smri'))
    numParameters = numParameters + 1;
    
    inputText(numParameters).promptString =  'How do you want to run Group ICA?';
    inputText(numParameters).answerString =  char('Serial', 'Parallel');
    inputText(numParameters).uiType = 'popup';
    inputText(numParameters).dataType = 'string';
    inputText(numParameters).tag = 'parallel_info';
    inputText(numParameters).enable = 'on';
    inputText(numParameters).value = 1;
    inputText(numParameters).flag = 'scalar';
    inputText(numParameters).help = struct('title', 'Run Parallel?', 'string', char('1. Serial - Group ICA will be run in serial mode. This option is preferred for analyzing smaller data-sets.', ...
        '2. Parallel - Group ICA will be run in parallel mode. Very useful option for analyzing large data. If you don''t have parallel computing toolbox parts of code will be run in sessions. Option is provided to enter number of sessions.'));
end


function matchInd = matchString(selectedStr, options, def)

if (ischar(selectedStr))
    matchInd = strmatch(lower(selectedStr), lower(cellstr(options)), 'exact');
else
    matchInd = selectedStr;
end

if (isempty(matchInd))
    matchInd = def;
end
