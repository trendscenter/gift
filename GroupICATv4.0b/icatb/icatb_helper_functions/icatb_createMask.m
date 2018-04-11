function mask_ind = icatb_createMask(files, fMRIHInfo, dataType, maskType, complexInfo)
%% Creates a mask based on the files given and mask type. If default mask is
% specified then first file of different subjects is loaded and indices greater than
% mean of all the voxels are taken into account. If mask is specified then
% only the indices in the mask are used.
%
% Input:
% 1. files - data structure with the field name
% 2. fMRIDIM - data dimensions
% 3. dataType - 'real' or 'complex'
% 4. complexInfo - structure containing file naming and complexType
%
% Output:
% 1. mask_ind - mask Indices


icatb_defaults;

global DEFAULT_MASK_OPTION;
global DEFAULT_MASK_SBM_MULTIPLIER;

dm_sbm_mult = DEFAULT_MASK_SBM_MULTIPLIER;
if (isempty(dm_sbm_mult))
    dm_sbm_mult = 0.01;
end

if ~exist('dataType', 'var')
    dataType = 'real';
end

if ~exist('maskType', 'var')
    maskType = 'default';
else
    if isempty(maskType)
        maskType = 'default';
    end
end

%% Default mask option
defaultMaskOption = DEFAULT_MASK_OPTION;

if isempty(defaultMaskOption)
    defaultMaskOption = 'first_file';
end

modalityType = icatb_get_modality;

%% Get the extension (checks only the first data set)
[pathstr, fName, extn] = fileparts(icatb_parseExtn(files(1).name(1, :)));

if (~strcmpi(extn, '.nii')  && ~strcmpi(extn, '.img') && ~strcmpi(extn, '.gz'))
    error('Mask should be specified in analyze or nifti format');
end

disp('Creating Mask');

% Initialise data count
countN = 0;

%% Mask type
if strcmpi(maskType, 'default')
    if strcmpi(modalityType, 'fmri')
        if strcmpi(defaultMaskOption, 'all_files')
            disp('Default mask includes voxels >= mean. Using all files of the subjects to create default mask ...');
        else
            disp('Default mask includes voxels >= mean. Using first file of each subject to create default mask ...');
        end
    elseif strcmpi(modalityType, 'smri')
        disp(['Using voxels >= ', num2str(100*dm_sbm_mult), '% of mean ...']);
    else
        disp('Using all EEG indices to create default mask ...');
    end
else
    disp(['Using file ', deblank(files(1).name(1, :)), ' to create mask ...']);
end


%% Default mask
if strcmpi(maskType, 'default')
    
    if strcmpi(modalityType, 'fmri')
        
        %% Loop over number of files
        for i = 1:length(files)
            
            %% Rename 4D nifti files
            tempF = icatb_rename_4d_file(files(i).name);
            
            %% Use first file or all files
            if ~strcmpi(defaultMaskOption, 'all_files')
                % use only the first time point
                tempF = deblank(tempF(1, :));
            end
            
            %% Include in-brain voxels
            for nn = 1:size(tempF, 1)
                countN = countN + 1;
                %% Load data
                data  = icatb_loadData(tempF, dataType, complexInfo, 'read', nn);
                if ~isreal(data)
                    data = abs(data);
                end
                
                % Add dummy baseline (Handle filtered data)
                data = data + 100;
                
                %% Use voxels that surpass or equal mean
                data = data(:);
                tempInd = (data >= mean(data));
                
                if (countN == 1)
                    nonZeroInd = tempInd;
                else
                    nonZeroInd = (nonZeroInd) & tempInd;
                end
                
                clear tempInd data;
            end
            %% End for including in-brain voxels
            
            clear tempF;
            
        end
        %% End for loop over number of files
        
    elseif strcmpi(modalityType, 'smri')
        
        % Rename 4D nifti files
        tempF = icatb_rename_4d_file(files(1).name);
        
        for nn = 1:size(tempF, 1)
            data = icatb_read_data(deblank(tempF(nn, :)));
            data = abs(data(:));
            tempInd = (data >= dm_sbm_mult*mean(data));
            if (nn == 1)
                nonZeroInd = tempInd;
            else
                nonZeroInd = (nonZeroInd) & tempInd;
            end
            
            clear tempInd data;
            
        end
    else
        nonZeroInd = ones(prod(fMRIHInfo.DIM), 1);
    end
    
else
    
    %% Calculate non-zero indices for selected mask
    dims = fMRIHInfo.V(1).dim(1:3); % Data dimensions
    data = icatb_loadData(deblank(files(1).name(1, :))); % Mask data
    size_data = size(data);
    if length(size_data) == 2
        size_data(3) = 1;
    end
    
    %% Check the dimensions of the mask w.r.t data
    if length(find(size_data == dims)) ~= length(dims)
        error('Error:MaskDim', 'Mask dimensions ([%s]) doesn''t match that of data dimensions ([%s])', ...
            num2str(size_data), num2str(dims));
    end
    
    %% Convert data to 2D matrix
    data = data(:);
    nonZeroInd = (data ~= 0);
    clear data;
    
end
%% End for checking default or selected mask

%% Mask out the zeros
mask_ind = find(nonZeroInd ~= 0);

if isempty(mask_ind)
    error('No voxels found. Error in creating Mask');
end

disp('Done Creating Mask');
