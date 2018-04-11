function [fileNames, zipFileName, files_to_zip] = icatb_saveICAData(outfile, icasig, A, mask_ind, numOfIC, HInfo, dataType, ...
    complexInfo, outputDir, zipFiles, deleteFiles)
% saves ICA data
%
% Input:
% 1. outfile - ouput file name
% 2. icasig - component images
% 3. A - time course
% 4. numOfIC - number of components
% 5. HInfo - header info
% 6. dataType - data type ('real' or 'complex')
% 7. complexInfo - complex Info
%
% Output:
%
% fileNames


if ~exist('mask_ind', 'var')
    error('mask indices must be passed');
end

if ~exist('dataType', 'var')
    dataType = 'real';
end

if ~exist('outputDir', 'var')
    outputDir = pwd;
end

if (~exist('complexInfo', 'var'))
    complexInfo = [];
end

%run defaults file
icatb_defaults;

global COMPONENT_NAMING;
global TIMECOURSE_NAMING;
global FUNCTIONAL_DATA_FILTER;
global ZIP_IMAGE_FILES;
global EEG_TOPOGRAPHY_NAMING;

if (~exist('zipFiles', 'var'))
    zipFiles = ZIP_IMAGE_FILES;
end

if (~exist('deleteFiles', 'var'))
    deleteFiles = 1;
end

[modalityType, dataTitle, compSetFields] = icatb_get_modality;

if ~strcmpi(modalityType, 'eeg')
    % get the image extension
    [pp, bb, fileExtn] = fileparts(FUNCTIONAL_DATA_FILTER);
else
    fileExtn = '.mat';
end

if size(A, 1) < size(A, 2)
    A = A';
end

if size(icasig, 1) ~= size(A, 2)
    icasig = icasig';
end


%determine output file names
lastUnderScore = icatb_findstr(outfile,'_');
lastUnderScore = lastUnderScore(end);
component_name = outfile(1:lastUnderScore);

if strcmpi(modalityType, 'eeg')
    eval([compSetFields{1}, ' = icasig;']);
    eval([compSetFields{2}, ' = A;']);
    clear icasig A;
    fileNames{1} = [component_name, fileExtn];
    icatb_save(fileNames{1}, compSetFields{:});
    zipFileName = {};
    files_to_zip = {};
    return;
end

%timecourse dimensions
tc_dim = [size(A, 1), numOfIC, 1];
%components dim, e.g. dimensions of image
c_dim = [HInfo.DIM(1),HInfo.DIM(2),HInfo.DIM(3)];

% Check the data type and if complex add R and I namings for the real part
% and imaginary part

%write timecourses in analyze format(1 file)
V = HInfo.V(1);
V.dim(1:3) = [tc_dim(1) tc_dim(2) 1];
V.dt(1) = 4;
V.n(1) = 1;

timecourse_name = strrep(component_name, COMPONENT_NAMING, TIMECOURSE_NAMING);

zipFileName = {};
if strcmpi(zipFiles, 'yes')
    % zip file name
    zipFileName = [component_name, '.zip'];
end


files_to_zip = {};

if isreal(A)
    % check if the time course is of double data type

    V.fname = [timecourse_name, fileExtn];

    % return zip files
    files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);

    V.fname = fullfile(outputDir, V.fname);

    % write the images
    icatb_write_vol(V, A);


elseif isa(A, 'complex_data')
    % check if it is complex data class

    % use file parts to separate time course name
    [pathstr, tName] = fileparts(timecourse_name);

    % File names of complex data
    P = icatb_get_complex_files_naming(timecourse_name, dataType, complexInfo, 'write');
    % get the file namings
    firstFileName = P.first; secondFileName = P.second;


    V.fname = firstFileName;

    % return zip files
    files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);

    V.fname = fullfile(outputDir, V.fname);

    % access the first field
    icatb_write_vol(V, getfield(A, 'firstField'));


    % phase image
    V.fname = secondFileName;

    % return zip files
    files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);

    V.fname = fullfile(outputDir, V.fname);

    % access the second field
    icatb_write_vol(V, getfield(A, 'secondField'));

    clear V1;

else

    error('Unknown data type');

end

% Initialise file names
fileNames = cell(1, numOfIC);

for i=1:numOfIC
    % return file index
    [fileIndex] = icatb_returnFileIndex(i);

    V = HInfo.V(1); V.dt(1) = 4;
    if strcmpi(fileExtn, '.img')
        % analyze format
        V.fname = [component_name, fileIndex, fileExtn];
        V.n(1) = 1;
    elseif strcmpi(fileExtn, '.nii')
        % handle nifti data
        % store the number in the component data
        V.n(1) = i;
        V.fname = [component_name, fileExtn];
    else
        error('Unknown image extensions');
    end

    fileNames{i} = V.fname;
    V.dim(1) = c_dim(1);
    V.dim(2) = c_dim(2);
    V.dim(3) = c_dim(3);

    data = zeros(prod(HInfo.DIM(1:3)), 1);
    if size(icasig, 2) == size(data, 1)
        data(mask_ind) = icasig(i, mask_ind);
    else
        data(mask_ind) = icasig(i, :);
    end
    % Reshape data
    data = reshape(data, [c_dim(1), c_dim(2), c_dim(3)]);
    %data = reshape(icasig(i, :), [c_dim(1), c_dim(2), c_dim(3)]);

    % handle complex data
    if isreal(data)
        if strcmpi(fileExtn, '.img')
            % return zip files
            files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);
        elseif strcmpi(fileExtn, '.nii')
            if i == numOfIC
                % return zip files
                files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);
            end
        end
        % full file name
        V.fname = fullfile(outputDir, V.fname);
        icatb_write_vol(V, data);
    elseif isa(data, 'complex_data')

        % File names of complex data
        P = icatb_get_complex_files_naming(V.fname, dataType, complexInfo, 'write');
        % get the file namings
        firstFileName = P.first; secondFileName = P.second;

        V.fname = firstFileName;

        if strcmpi(fileExtn, '.img')
            % return zip files
            files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);
        elseif strcmpi(fileExtn, '.nii')
            if i == numOfIC
                % return zip files
                files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);
            end
        end

        V.fname = fullfile(outputDir, V.fname);
        % write image
        icatb_write_vol(V, getfield(data, 'firstField'));

        V.fname = secondFileName;

        if strcmpi(fileExtn, '.img')
            % return zip files
            files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);
        elseif strcmpi(fileExtn, '.nii')
            if i == numOfIC
                % return zip files
                files_to_zip = returnZipFiles(V.fname, fileExtn, files_to_zip, zipFiles);
            end
        end

        V.fname = fullfile(outputDir, V.fname);
        % write image
        icatb_write_vol(V, getfield(data, 'secondField'));

    else

        error('Unknown data type.');

    end

end
% end for all components

% zip files
if (~isempty(zipFileName))
    [p, fN, extn] = fileparts(zipFileName);
    outputDir = fullfile(outputDir, p);
    zipFileName2 = [fN, extn];
    files_to_zip2 = regexprep(files_to_zip, ['.*\', filesep], '');
    icatb_zip(zipFileName2, files_to_zip2, outputDir);
    if (deleteFiles)
        % delete the files
        icatb_delete_file_pattern(char(files_to_zip2), outputDir);
    end
end
% end for checking

function files_to_zip = returnZipFiles(newFile, fileExtn, files_to_zip, zipFiles)

if ~exist('files_to_zip', 'var')
    files_to_zip = {};
end

if strcmpi(zipFiles, 'yes')
    countZip = length(files_to_zip) + 1;
    % for analyze zip header files also
    if strcmpi(fileExtn, '.img')
        files_to_zip{countZip} = [newFile];
        countZip = length(files_to_zip) + 1;
        files_to_zip{countZip} = [newFile(1:end-3), 'hdr'];
    else
        files_to_zip{countZip} = [newFile];
    end

end
