function [data, HInfo, XYZ] = icatb_read_data(P, file_numbers, mask_ind, precisionType, varToLoad)
%% Function to read the data using spm_vol and spm_read_vols functions.
%
% Inputs:
% 1. P - Character array of file names (Use full path)
% 2. file_numbers - File number/numbers to include
% 3. mask_ind - Mask indices
% 4. precisionType - Double or Single precision
%
% Outputs:
% 1. data - 4D double array (x, y, z, t) and if mask is specified for fMRI 2D double array (Voxels, t) is
% returned.
% 2. HInfo - Header information
% 3. XYZ - Real world coords
%

%% Use defaults
if (iscell(P))
    P = char(P);
end

if (~exist('file_numbers', 'var'))
    file_numbers = [];
end

if (~exist('mask_ind', 'var'))
    mask_ind = [];
end

if (~isempty(mask_ind) && ischar(mask_ind))
    mV = icatb_read_data(mask_ind);
    clear mask_ind;
    mask_ind = find(abs(mV) > eps);
    clear mV;
end

if (~exist('precisionType', 'var'))
    precisionType = 'double';
end

if (~exist('varToLoad', 'var'))
    varToLoad = '';
end


if (~isempty(file_numbers))
    P = icatb_rename_4d_file(P);
end

%% Parse extension
first_file = deblank(P(1, :));
first_file = icatb_parseExtn(first_file);
[pp, fileN, extn] = fileparts(first_file);


%% For EEG data
if strcmpi(extn, '.mat')
    
    XYZ = []; % Return real world coordinates as empty
    
    if ~isempty(file_numbers)
        P = P(file_numbers, :);
    end
    
    %% Loop over files
    for nP = 1:size(P, 1)
        % Load data
        if (isempty(varToLoad))
            dataInfo = load(deblank(P(nP, :)));
        else
            dataInfo = load(deblank(P(nP, :)), varToLoad);
        end
        fieldN = fieldnames(dataInfo);
        temp = getfield(dataInfo, fieldN{1});
        if ~isreal(temp)
            error(['Check file ', deblank(P(nP, :)), ' as it is not a numeric data type']);
        end
        
        if (length(size(temp)) == 2)
            temp = reshape(temp, size(temp, 1), 1, size(temp, 2));
        end
        
        % temp = squeeze(temp);
        
        % Initialise data
        if nP == 1
            if (~isempty(mask_ind))
                if strcmpi(precisionType, 'double')
                    data = zeros(length(mask_ind), size(temp, 3), size(P, 1));
                else
                    data = zeros(length(mask_ind), size(temp, 3), size(P, 1), 'single');
                end
            else
                if strcmpi(precisionType, 'double')
                    data = zeros(size(temp, 1), size(temp, 2), 1, size(temp, 3), size(P, 1));
                else
                    data = zeros(size(temp, 1), size(temp, 2), 1, size(temp, 3), size(P, 1), 'single');
                end
            end
        end
        % End for initialising data
        
        if strcmpi(precisionType, 'single')
            temp = single(temp);
        end
        
        if (nP == 1)
            HInfo.DIM = [size(temp, 1), size(temp, 2), 1];
            HInfo.V.n(1) = 1;
            HInfo.V.dim = HInfo.DIM;
        end
        
        if (~isempty(mask_ind))
            temp = reshape(temp, prod(HInfo.DIM), size(temp, length(size(temp))));
            data(:, :, nP) = temp(mask_ind, :);
        else
            data(:, :, 1, :, nP) = temp;
        end
        %clear temp
    end
    %% End loop over files
    
    %     if (~isempty(mask_ind))
    %         data = reshape(data, prod(HInfo.DIM), size(data, 4), size(data, 5));
    %         data = squeeze(data(mask_ind, :, :));
    %     end
    
    return;
end


warning('off', 'MATLAB:dispatcher:ShadowedMEXExtension');

if (~strcmpi(extn, '.gz'))
    
    if ~isempty(file_numbers)
        try
            %% Try getting the volume of the given file
            V = icatb_spm_vol(deblank(P(file_numbers, :)));
        catch
            %% Get full volume (Handle this for the files that were selected
            % using the old version)
            V = icatb_spm_vol(P);
            if (max(file_numbers) > length(V))
                error('Error:Timepoints', 'Maximum of file numbers (%d) specified exceeds the no. of time points (%d)\n', ...
                    max(file_numbers), length(V));
            end
            % Truncate the volume
            V = V(file_numbers);
        end
    else
        V = icatb_spm_vol(P);
    end
    % end for checking the existence of file number variable
    
    if (~isempty(mask_ind))
        if strcmpi(precisionType, 'double')
            %% Use less memory if mask is specified
            data = zeros(length(mask_ind), length(V));
        else
            data = zeros(length(mask_ind), length(V), 'single');
        end
    else
        if strcmpi(precisionType, 'double')
            %% 4D double array
            data = zeros([V(1).dim(1:3), length(V)]);
        else
            data = zeros(V(1).dim(1), V(1).dim(2), V(1).dim(3), length(V), 'single');
        end
    end
    
    if (~isempty(mask_ind))
        voxel_coords = icatb_get_voxel_coords(V(1).dim(1:3));
        voxel_coords = voxel_coords(:, mask_ind);
    end
    
    
    %% Loop over time points
    for nVol = 1:length(V)
        
        %% Check dimensions
        if (nVol == 1)
            oldDim = V(nVol).dim(1:3);
        else
            newDim = V(nVol).dim(1:3);
            if (length(find((newDim == oldDim) == 1)) ~= 3)
                error('Error:DIM', 'Please check file %s as this doesn''t equal the dimension [%d, %d, %d]\n', V(nVol).fname, ...
                    newDim(1), newDim(2), newDim(3));
            end
        end
        
        if ((nVol == 1) && (nargout == 3))
            %% Read data time point by time point and apply mask
            [tempData, XYZ] = icatb_spm_read_vols(V(nVol));
            if (~isempty(mask_ind))
                tempData = tempData(mask_ind);
            end
        else
            if (~isempty(mask_ind))
                tempData = icatb_spm_get_data(V(nVol), voxel_coords, 0);
            else
                tempData = icatb_spm_read_vols(V(nVol));
            end
        end
        
        tempData(isfinite(tempData) == 0) = 0;
        
        if strcmpi(precisionType, 'single')
            tempData = single(tempData);
        end
        
        if (~isempty(mask_ind))
            data(:, nVol) = tempData(:);
        else
            data(:, :, :, nVol) = tempData;
        end
        
        clear tempData;
        
    end
    %% End of loop over time points
    
    %% Get the voxel size
    VOX = double(V(1).private.hdr.pixdim(2:4)); VOX = abs(VOX);
    
    
else
    
    [data, hdr, V] = loadNIIGZ(P, file_numbers, mask_ind);
    
    if (nargout == 3)
        XYZ = icatb_get_voxel_coords(V(1).dim(1:3));
        RCP = XYZ;
        RCP(4, :) = 1;
        XYZ = (V(1).mat(1:3, :)*RCP);
    end
    
    VOX = abs(double(hdr.pixdim(2:4)));
    
end

if ((nargout == 3) && ~isempty(mask_ind))
    XYZ = XYZ(:, mask_ind);
end

if (nargout >= 2)
    %% Save Header info in structure
    HInfo = struct('DIM', (V(1).dim(1:3)), 'V', V, 'VOX', VOX);
end

clear V;


function [data, hdr, V] = loadNIIGZ(P, file_numbers, mask_ind)

P = icatb_rename_4d_file(P);

if (isempty(file_numbers))
    file_numbers = (1:size(P, 1));
end

P = P(file_numbers, :);
P = cellstr(P);

fNames = regexprep(P, '(.*),(.*)','$1');
fnums = regexprep(P, '(.*),(.*)','$2');

[dd, inds] = unique(fNames);
inds = sort(inds);
ufNames = fNames(inds);
dat = cell(length(inds), 1);
for nI = 1:length(inds)
    chk = strcmp(fNames, ufNames{nI});
    tmp_inds = find(chk ~= 0);
    nums = str2num(char(fnums(tmp_inds(1):tmp_inds(end))));
    [dat{nI}, hdr, tmpV] = icatb_read_gzip_nii(ufNames{nI}, 'timepoints', nums, 'mask', mask_ind);
    if (nI == 1)
        V = tmpV;
    end
end

data = cat(length(size(dat{1})), dat{:});


