function sesInfo = icatb_update_mask(sesInfo)
% Create mask from the data

icatb_defaults;
global REMOVE_CONSTANT_VOXELS;

removeConstVoxels = 0;
if (~isempty(REMOVE_CONSTANT_VOXELS))
    removeConstVoxels = REMOVE_CONSTANT_VOXELS;
end

[sesInfo, complexInfo] = icatb_name_complex_images(sesInfo, 'read');

% Get modality type
[modalityType, dataTitle] = icatb_get_modality;

default_mask_opts = icatb_default_mask_opts;

parallelMode = 'serial';
num_workers = 4;
try
    parallelMode = sesInfo.userInput.parallel_info.mode;
catch
end

try
    num_workers = sesInfo.userInput.parallel_info.num_workers;
catch
end

toolboxNames = ver;
parallelCluster = ~isempty(find(strcmpi(cellstr(char(toolboxNames.Name)), 'parallel computing toolbox') == 1));
runParallel = 0;
if (strcmpi(parallelMode, 'parallel') && parallelCluster)
    runParallel = 1;
end

if (runParallel)
    try
        matlabpool('open', num_workers);
    catch
    end
end


if ~strcmpi(modalityType, 'eeg')
    % Handle fmri data here
    
    % return volume and header information
    [V, HInfo] = icatb_returnHInfo(deblank(sesInfo.userInput.files(1).name(1, :)));
    V = V(1);
    
    %get dimensions of data
    xdim = HInfo.DIM(1); ydim = HInfo.DIM(2); zdim = HInfo.DIM(3);
    
    % assign HInfo to sesInfo.userInput
    sesInfo.userInput.HInfo = HInfo;
    
    if ~isfield(sesInfo.userInput, 'maskFile')
        sesInfo.userInput.maskFile = [];
    end
    
    % Calculate mask
    if(isempty(sesInfo.userInput.maskFile))
        flagMask = 'default';
        if (~runParallel)
            mask_ind = icatb_createMask(sesInfo.userInput.files, HInfo, ...
                sesInfo.userInput.dataType, flagMask, complexInfo);
        else
            mask_ind = icatb_parCreateMask(sesInfo.userInput.files, HInfo, ...
                sesInfo.userInput.dataType, flagMask, complexInfo);
        end
    else
        flagMask = 'user_specified';
        files(1).name = sesInfo.userInput.maskFile;
        mask_ind = icatb_createMask(files(1), HInfo, 'real', flagMask, complexInfo);
        clear files;
    end
    
    mask = zeros(xdim, ydim, zdim);
    
    if (strcmpi(modalityType, 'fmri') && removeConstVoxels)
        nonZeroInd = mask;
        nonZeroInd(mask_ind) = 1;
        %XYZ = icatb_get_voxel_coords(HInfo.DIM(1:3));
        mask_ind = remove_constant_voxels(sesInfo.userInput.files, nonZeroInd(:), flagMask);
    end
    
    % Write mask image
    mask(mask_ind) = 1;
    V.fname = fullfile(sesInfo.userInput.pwd, [sesInfo.userInput.prefix, 'Mask.img']);
    V.n(1) = 1;
    icatb_write_vol(V, mask);
    
    
else
    % Handle EEG data here
    [origData, HInfo] = icatb_loadData(deblank(sesInfo.userInput.files(1).name(1, :)));
    sesInfo.userInput.HInfo = HInfo;
    % Data dimensions
    dims = [size(origData, 1), size(origData, 2), size(origData, 4)];
    
    if (isempty(sesInfo.userInput.maskFile))
        
        % Use all indices
        temp = ones(prod(dims(1:2)), 1);
    else
        
        files(1).name = sesInfo.userInput.maskFile;
        
        temp = icatb_loadData(files(1).name);
        
        if length(find(size(temp) == dims(1:2))) ~= 2
            error(['Mask ', files(1).name, ' doesn''t match data dimensions ']);
        end
        
    end
    % Mask indices
    mask_ind = find(temp(:) ~= 0);
    
end
% End for checking data type

% Update sesInfo variable
sesInfo.userInput.mask_ind = mask_ind;

if (exist('default_mask_opts', 'var'))
    sesInfo.userInput.default_mask_opts = default_mask_opts;
end

function mask_ind = remove_constant_voxels(files, nonZeroInd, flagMask)
%% Remove constant voxels
%

errMsg = 'No non-constant voxels found. Error in creating Mask';
fprintf('\n');
disp('Removing constant voxels from all subjects ...');
mask_inds = find(nonZeroInd == 1);
for i = 1:length(files)
    
    dat = read_data(files(i).name, mask_inds);
    % Add dummy baseline
    dat = dat + 100;
    tmpM = any(diff(dat, 1));
    if (i == 1)
        chk = tmpM;
    else
        chk = chk & tmpM;
    end
    
    clear tmpM;
    
    if (isempty(find(chk == 1)))
        error(errMsg);
    end
    
    if (~strcmpi(flagMask, 'default'))
        try
            tmpM = (prod(double(dat >= repmat(0.01*mean(dat, 2), 1, size(dat, 2)))) ~= 0);
        catch
            % Remove voxels at the edges
            for nD = 1:size(dat, 1)
                tmp = (dat(nD, :) >=  0.01*mean(dat(nD, :)));
                if (nD == 1)
                    tmpM = tmp;
                else
                    tmpM = tmpM & tmp;
                end
            end
        end
        chk = chk & tmpM;
    end
    
    clear dat tmpM;
end

nonZeroInd(nonZeroInd == 1) = chk(:);
mask_ind = find(nonZeroInd ~= 0);
if (isempty(mask_ind))
    error(errMsg);
end

disp('Done');
fprintf('\n');



function Y = read_data(files, mask_ind)
%% Read data
%

Y = icatb_read_data(files, [], mask_ind);
Y = Y';
%Y = icatb_spm_get_data(icatb_spm_vol(files), XYZ, 0);
%Y(isfinite(Y) == 0) = 0;
