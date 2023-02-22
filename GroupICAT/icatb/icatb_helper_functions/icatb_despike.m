function icatb_despike(files, TR, maskFile, outputDir)
%% Despike fmri timeseries
%
% Inputs:
% 1. files - Cell array of file names
% 2. TR - Experimental TR in seconds
% 3. Maskfile - Mask file name.
% 4. outputDir - Output directory to place results
%
%
% Output files are written to directory outputDir with d_sub* prefix.
%

if (~exist('files', 'var'))
    files = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'filter', '*nii', 'title', 'Select image files for despiking ...');
end

drawnow;

if (isempty(files))
    error('Files are not selected.');
end

if (~exist('TR', 'var'))
    TR = icatb_inputdlg2('Enter TR in seconds', 'TR', 1, {''});
    if (isempty(TR))
        error('TR is not selected');
    end
    TR = str2num(TR{1});
end

if (isempty(TR))
    error('TR specified is not a valid number');
end


if (~exist('maskFile', 'var'))
    chkMask = icatb_questionDialog('title', 'Mask?', 'textbody', 'Do You Want To Select Mask?');
    maskFile = '';
    if (chkMask)
        maskFile = icatb_selectEntry('typeSelection', 'single', 'typeEntity', 'file', 'filter', '*img;*nii', 'filetype', 'image', 'title', 'Select mask. Dimensions must match the data.');
        if (isempty(maskFile))
            error('mask is not selected');
        end
    end
end

if (~exist('outputDir', 'var'))
    outputDir = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'directory', 'title', 'Select output directory ...');
end

if (isempty(outputDir))
    error('Directory is not selected.');
end

files = cellstr(files);

files = getFullPaths(files);

disp('Despiking data. This process might take a while ...');

drawnow;

for nF = 1:length(files)
    
    cF = files{nF};
    
    [pp, fn] = fileparts(deblank(cF(1, :)));
    
    disp(['Processing Subject ', cF(1, :), ' ...']);
    
    if (isempty(maskFile))
        [data, HInfo] = icatb_loadData(cF);
        dims = [size(data, 1), size(data, 2), size(data, 3), size(data, 4)];
        data = reshape(data, prod(dims(1:3)), dims(4))';
    else
        mask = icatb_loadData(maskFile);
        mask = squeeze(mask(:, :, :,1));
        mask_inds = find(abs(mask) > eps);
        [data, HInfo] = icatb_read_data(cF, [], mask_inds);
        dims = [HInfo.DIM(1:3), size(data, 2)];
        data = data';
    end
    
    [data, frames_replaced] = icatb_despike_tc(data, TR);
    if (isempty(maskFile))
        data = reshape(data', dims);
    else
        dataN = zeros(prod(dims(1:3)), dims(4));
        dataN(mask_inds, :) = data';
        data = reshape(dataN, dims);
        clear dataN;
    end
    
    outputFileName = fullfile(outputDir, ['d_sub_', icatb_returnFileIndex(nF), '_', fn, '.nii']);
    
    disp(['Writing out despiked data in file ', outputFileName, ' ...']);
    
    icatb_write_nifti_data(outputFileName, HInfo.V, data, 'despiked data');
    
    % frames replaced
    frameOutFname = fullfile(outputDir, ['d_sub_', icatb_returnFileIndex(nF), '_', fn, '_frames_info.mat']);
    save(frameOutFname, 'frames_replaced');
    
    drawnow;
    
end

fprintf('\nDone\n');


function files = getFullPaths(files)
%% Get full file paths

for i = 1:length(files)
    
    fn = files{i};
    if (size(fn, 1) == 1)
        [p, fN, extn] = fileparts(fn);
        f = dir(fn);
        f = f([f.isdir] == 0);
        if (isempty(f))
            error(['No files found. Please check the file pattern of ', fn]);
        end
        if ~strcmp(p(end), filesep)
            p = [p, filesep];
        end
        files{i} = [repmat(p, length(f), 1), char(f.name)];
    end
end
