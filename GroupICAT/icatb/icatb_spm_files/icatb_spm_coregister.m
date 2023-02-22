function resliced_files = icatb_spm_coregister(reference, source, other_pattern)

defaults = spm_get_defaults;

orig_dir = pwd;
%cd(directory);
%progFile = fullfile(orig_dir,'cs_progress.txt');

%cs_log( ['Beginning cs_coregister for ', pwd], progFile );

% Source image
%[coregSourceImage] = check_file_path(source, 'source');
coregSourceImage = source;

otherImages = other_pattern;

% if ~isempty(other_pattern)
%     % Other images used
%     otherImages = cs_list_files(pwd, other_pattern, 'fullpath');
% else
%     otherImages = '';
% end
% 
% if ~isempty(otherImages)
%     otherImages = icatb_rename_4d_file(otherImages);
% else
%     otherImages = '';
% end


% Coregister step
%if (csprefs.run_coreg)

% Reference image
%[coregRefImage] = check_file_path(reference, 'reference');
coregRefImage = reference;

% Estimate options pulled from spm_defaults file
eoptions = defaults.coreg.estimate;

% The above code is from spm_config_coreg function
x  = spm_coreg(coregRefImage, coregSourceImage, eoptions);

M  = inv(spm_matrix(x));

PO = strvcat(strvcat(coregSourceImage), strvcat(otherImages));

MM = zeros(4,4,size(PO,1));

for j=1:size(PO,1),
    MM(:,:,j) = spm_get_space(deblank(PO(j,:)));
end;

for j=1:size(PO,1),
    spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
end;

%cs_log(['spm_coreg completed for ', pwd], progFile);

%end

% Reslice step
%if (csprefs.run_reslice)

% Reference image
%[refImage] = check_file_path(reference, 'reference');
refImage = coregRefImage;

P = strvcat(strvcat(refImage), strvcat(strvcat(coregSourceImage), strvcat(otherImages)));

% Reslice options used from spm_defaults
flags.mask   = defaults.coreg.write.mask;
flags.mean   = 0;
flags.interp = defaults.coreg.write.interp;
flags.which  = 1;
flags.wrap   = defaults.coreg.write.wrap;

spm_reslice(P, flags);

%cs_log(['spm_reslice completed for ', pwd], progFile);

%end

resliced_files = cell(size(otherImages, 1), 1);
for nF = 1:size(otherImages, 1)
    [p, fn, extn] = fileparts(deblank(otherImages(nF,:)));
    resliced_files{nF} = fullfile(p, ['r', fn, extn]);
end

%cs_log( ['Ending cs_coregister for ', pwd], progFile );

cd(orig_dir);



function [files] = check_file_path(filePattern, optional)
% Form files from a file pattern and select only the first file

if exist('optional', 'var')
    optional = '';
end

% check if full file path is specified or it is a file on Matlab path
if (exist(filePattern, 'file') == 2)
    files = filePattern;
else
    files = cs_list_files(pwd, filePattern, 'fullpath');
end

if isempty(files)
    error(['Please check the file pattern: ', filePattern, ' for ', optional, ' image']);
end

% Add a number at the end for 4D Nifti files
files = icatb_rename_4d_file(files);

if size(files, 1) > 1
    files = deblank(files(1, :));
end

function files = cs_list_files(inputDir, filePattern, optional)
%% List files using the file pattern
%
% Input:
% 1. inputDir - Input directory
% 2. filePattern - File pattern
% 3. optional - Optional variable returns relative or full file path
%
% Ouput:
% files - files with the specified file pattern


% Use regular expression for files
global FILE_USEREGEXP;

files = '';

if ~exist('inputDir', 'var')
    error('Input directory is not specified for listing file pattern');
end

if ~exist('filePattern', 'var')
    filePattern = '*';
end

if ~exist('optional', 'var')
    optional = 'relative';
end

%% Don't use regular expression if empty
if isempty(FILE_USEREGEXP)
    FILE_USEREGEXP = 0;
end

%% Change directory to get inputdirectory name
oldDir = pwd;

cd(inputDir);

%% Handle nested file pattern
if (~FILE_USEREGEXP)
    [rel_dir, filePattern, extn] = fileparts(filePattern);
    filePattern = [filePattern, extn];
    if (~isempty(rel_dir))
        try
            cd(rel_dir);
        catch
            if (ispc)
                rel_dir = dir(rel_dir);
            else
                rel_dir = ls('-d', rel_dir);
            end
            cd(deblank(rel_dir));
        end
    end
end

inputDir = pwd;

cd(oldDir);
% End for changing directory

% List files using dir function
if ~FILE_USEREGEXP
    d = dir(fullfile(inputDir, filePattern));
else
    d = dir(inputDir);
end

isDirs = [d.isdir];

checkIndices = find(isDirs ~= 1);

if ~isempty(checkIndices)
    
    files = str2mat(d(checkIndices).name);
    
    if FILE_USEREGEXP
        % Use regular expression
        files = cellstr(str2mat(d(checkIndices).name));
        if ~ispc
            inds = regexp(files, filePattern);
        else
            inds = regexpi(files, filePattern);
        end
        inds = cs_good_cells(inds);
        if ~isempty(find(inds))
            files = str2mat(files(inds));
        else
            clear files;
            files = '';
            return;
        end
    end
    
    if ~strcmpi(optional, 'relative')
        
        if ~strcmp(inputDir(end), filesep)
            inputDir = [inputDir, filesep];
        end
        
        % Return full file path
        files = strcat(inputDir, files);
        
    end
    
end
