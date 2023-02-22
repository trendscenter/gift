function [icasig, icaTimecourse, structuralImage, coords, HInfo, text_left_right, classComp] = icatb_loadICAData(varargin)
% load ICA time course
%
% Inputs must be in pairs:
%
% 1. structFile - structural file
% 2. compFiles - component files
% 3. slicePlane - anatomical plane
% 4. sliceRange - slices in mm
% 5. comp_numbers - component numbers
% 6. returnValue - image values
% 7. convertToZ - convert to zscores
% 8. threshValue - threshold value
%
% Output:
% icasig - component image
% icaTimecourse - ica Time course
% structural Image - structural image
% coords - Real world coordinates
% HInfo - header info
% text_left_right - Text for plotting right and left
% classComp - Check if the component images are of 3d analyze or nifti
% class

icatb_defaults;
global IMAGE_VALUES;
global CONVERT_Z;
global THRESHOLD_VALUE;
global IMAGES_PER_FIGURE;
global ANATOMICAL_PLANE;
global ZIP_IMAGE_FILES;
global COMPONENT_NAMING;
global TIMECOURSE_NAMING;
global FLIP_ANALYZE_IMAGES;

% image values
imageValues = IMAGE_VALUES;

% imageReshape
imReshape = 'yes';

openDialog = 'yes';

dataType = 'real';

complexInfo = [];

zipFileName = [];

files_in_zip = {};

mask_file = [];

flag_delete = 'yes';

% get the image value
switch lower(imageValues)
    case 'positive'
        returnValue = 2;
    case 'absolute value'
        returnValue = 3;
    case 'negative'
        returnValue = 4;
    otherwise
        returnValue = 1;
end

% get the threshold value
threshValue = eval(THRESHOLD_VALUE);

% anatomical plane
slicePlane = lower(ANATOMICAL_PLANE);

% convert to z scores
if strcmpi(CONVERT_Z, 'yes')
    convertToZ = 1;
else
    convertToZ = 0;
end

flip_analyze_images = [];
%display_method = 'component explorer';

% get the required parameters
for ii = 1:2:nargin

    if strcmpi(varargin{ii}, 'structfile')
        % structural file
        structFile = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'compfiles')
        % component files
        compFiles = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'sliceplane')
        % slice plane
        slicePlane = lower(varargin{ii + 1});
    elseif strcmpi(varargin{ii}, 'slicerange')
        % slice range
        sliceRange = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'comp_numbers')
        % component numbers
        comp_numbers = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'returnvalue')
        % return value
        returnValue = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'converttoz')
        % convert to z scores
        convertToZ = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'threshvalue')
        % threshold value
        threshValue = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'imreshape')
        % reshape image to 3D
        imReshape = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'open_dialog')
        % open help dialog
        openDialog = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'complexinfo')
        % complex info
        complexInfo = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'datatype')
        % data type
        dataType = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'zipfile')
        zipFileName = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'files_in_zip_file')
        files_in_zip = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'mask_file')
        mask_file = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'flag_delete')
        flag_delete = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'flip_analyze_images')
        flip_analyze_images = varargin{ii + 1};
    end
    % end for checking parameters

end
% end for loop

if ~exist('structFile', 'var')
    error('Structural file is missing');
end

if ~exist('compFiles', 'var')
    error('Component files are missing');
end

% new_compFiles = get_complex_files_naming(compFiles, dataType, complexInfo, 'write');
% if isstruct(new_compFiles)
%     new_compFiles = str2mat(new_compFiles.first);
% end
% numComp = icatb_get_countTimePoints(new_compFiles);
numComp = length(comp_numbers);

% % get the default component numbers
% if ~exist('comp_numbers', 'var')
%     comp_numbers = [1:1:numComp];
% end

if ~exist('sliceRange', 'var')
    % get the slices in mm
    slicePara = icatb_get_slice_def(icatb_get_vol_nifti(structFile), slicePlane);
    sliceRange = slicePara.slices;
    clear slicePara;
end

% Open help dialog
if strcmpi(openDialog, 'yes')
    % resize the image and return header Info
    helpHandle = helpdlg('Interpolating component images. Please wait ...', 'Interpolating components');
end
% end for comparing help dialog


[pathstr_comp, c_name, compExtn] = fileparts(deblank(compFiles(1, :)));

% Unzip files
if ~isempty(zipFileName)
    icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), pathstr_comp);
end

%%%%%%%% Display warning message %%%%%%%%%%
compV = icatb_get_vol_nifti(deblank(compFiles(1, :)));
classComp = 'analyze';
if strcmpi(class(compV(1).private), 'icatb_nifti')
    classComp = 'nifti';
end

[status] = checkMatFile(deblank(compFiles(1, :)));

if ~isempty(flip_analyze_images)
    if strcmpi(classComp, 'analyze') && (status == 0)
        if (flip_analyze_images ~= FLIP_ANALYZE_IMAGES)
            warning(['Flip parameter (', num2str(FLIP_ANALYZE_IMAGES), ...
                ') specified in icatb_defaults is different from the one specified during the analysis (', ...
                num2str(flip_analyze_images), ')']);
        end
    end
end
%%%%%%%% End for displaying warning message %%%%%%%%%%

%%%%%%%%% Resize the image to structural dimensions and apply the
%%%%%%%%% necessary transformations %%%%%%
% resize the image and return header Info
[structuralImage, icasig, coords, HInfo, text_left_right] = icatb_returnResizedImage(structFile, compFiles, ...
    slicePlane, sliceRange, dataType, complexInfo, comp_numbers, mask_file);

clear parameters;

%%%%%%%%% Apply structural image parameters %%%%%%%%
structHInfo = HInfo;
structDIM = HInfo.DIM;

icasig = reshape(icasig, size(icasig, 1), prod(structDIM));

% apply display parameters
[icasig] = icatb_applyDispParameters(icasig, convertToZ, returnValue, threshValue, ...
    structDIM, HInfo);

% reshape image to 3D
if strcmpi(imReshape, 'yes')
    icasig = reshape(icasig, [size(icasig,1), HInfo.DIM(1), HInfo.DIM(2), HInfo.DIM(3)]);
end
icaDIM = HInfo.DIM;

% load time course
icaTimecourse = icatb_loadICATimeCourse(compFiles, dataType, complexInfo, comp_numbers);

if strcmpi(flag_delete, 'yes')
    if ~isempty(zipFileName)
        files_in_zip = str2mat(regexprep(files_in_zip, ['.*\', filesep], ''));
        %fileN = icatb_fullFile('files', files_in_zip, 'directory', pathstr_comp);
        icatb_delete_file_pattern(files_in_zip, pathstr_comp);
    end
end

try
    delete(helpHandle);
    clear helpHandle;
catch
end


function [status] = checkMatFile(file_name)
% Check MAT file name and return 1 if it exists

status = 0;
file_name = icatb_parseExtn(file_name);
[pp, fName, extn] = fileparts(file_name);

% MAT file name
if strcmp(extn, '.img')
    matFileName = fullfile(pp, [fName, '.mat']);
else
    matFileName = fullfile(pp, [fName, '.MAT']);
end

status = exist(matFileName, 'file');

if status > 1
    status = 1;
end