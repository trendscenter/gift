function icatb_talairach(files, threshold, convertToZ, sample_size)
%% Write talariach tables using Dae Il Kim's script and talariach daemon.
% 1. file - Input image file/files
% 2. threshold - Default is 3.5.
% 3. convertToZ - Options are 0 and 1.
% 4. sample_size - Distance between contiguous voxels in mm.
%

% Defaults
icatb_defaults;
global TALAIRACH_DIST;
global TALAIRACH_THRESHOLD;

isGUI = 0;
if ~exist('files', 'var')
    isGUI = 1;
    files = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'fileType', 'image', 'title', 'Select image/images ...', 'filter', ...
        '*.img;*.nii');
end

drawnow;

if isempty(files)
    error('Input image file/files is/are not selected');
end

if ~exist('threshold', 'var')
    threshold = TALAIRACH_THRESHOLD;
end

if (~exist('sample_size', 'var'))
    sample_size = TALAIRACH_DIST;
end

if (~exist('convertToZ', 'var'))
    convertToZ = 0;
end

if (isGUI)
    [threshold, sample_size, convertToZ] = getParams(threshold, sample_size);
end

drawnow;

try
    setJVM;
catch
end

tal_prog_file = fullfile(fileparts(which('gift.m')), 'icatb_talairach_scripts', 'talairach', 'talairach.jar');

% Run talairach daemon without internet (talairach.jar file updated)
commandStr = ['java -classpath "', tal_prog_file, '" org.talairach.AtlasServer 1600&'];

[status, result] = system(commandStr);

if (status == 1)
    error(result);
end

drawnow;


files = cellstr(files);

files = formFullPaths(files);

threshold = abs(threshold);

disp(['Selected threshold is ', num2str(threshold)]);

%% Loop over files
for nFiles = 1:length(files)
    
    file = files{nFiles};
    
    fprintf('\n');
    
    msgStr = ['Loading ', file, ' ...'];
    disp(msgStr);
    
    if (convertToZ)
        disp('Converting to z-scores ...');
    end
    
    V = icatb_spm_vol(file);
    
    %% Loop over no. of volumes
    for nV = 1:length(V)
        
        [outputDir, file_name] = fileparts(V(nV).fname);
        
        if ((length(V) > 1) || (V(nV).n(1) > 1))
            file_name = [file_name, '_', icatb_returnFileIndex(V(nV).n(1))];
        end
        
        data = icatb_spm_read_vols(V(nV));
        data(isnan(data)) = 0;
        
        %% Convert to z-scores
        if (convertToZ)
            tmp = data;
            tmp = tmp(:);
            mask_indices = (tmp ~= 0);
            tmp(mask_indices) = detrend(tmp(mask_indices), 0)./std(tmp(mask_indices));
            data = reshape(tmp, V(nV).dim(1:3));
            clear tmp
        end
        
        
        % Voxel dimensions %
        xdim = V(nV).dim(1); ydim = V(nV).dim(2); zdim = V(nV).dim(3);
        calcOrigin = (V(nV).mat\[0 0 0 1]')';
        origin = calcOrigin(1:3);
        voxel1 = V(nV).mat(1, 1);
        voxel2 = abs(V(nV).mat(2, 2));
        voxel3 = abs(V(nV).mat(3, 3));
        
        %% Report Left Right in Neurological convention
        volStr = 'volume (cc)';
        valStr = 'random effects: Max Value (x, y, z)';
        
        drawnow;
        
        fprintf('\n');
        
        disp('Writing talaraich coordinates for positive values ...');
        fprintf('\n');
        write_tal('pos', outputDir, file_name, data, sample_size, threshold, xdim, ydim, zdim, voxel1, voxel2, voxel3, origin, volStr, valStr);
        fprintf('\n');
        disp('Done writing talaraich coordinates for positive values');
        fprintf('\n');
        
        disp('Writing talaraich coordinates for negative values ...');
        fprintf('\n');
        write_tal('neg', outputDir, file_name, data, sample_size, threshold, xdim, ydim, zdim, voxel1, voxel2, voxel3, origin, volStr, valStr);
        fprintf('\n');
        disp('Done writing talaraich coordinates for negative values');
        fprintf('\n');
        
    end
    %% End of loop over no. of volumes
    
end
%% End of loop over files


function write_tal(flagCheck, outputDir, file_name, data, sample_size, threshold, xdim, ydim, zdim, voxel1, voxel2, voxel3, origin, volStr, valStr)
% Write talairach coords

TALAIRACH_FILE = 'talairach';
TALAIRACH_EXCEL = 'talairach_final';

if ~strcmpi(flagCheck, 'pos')
    data = -data;
end

% voxel size dependent on sample_size
voxsize = (sample_size^3)/1000;

% Location of the Talairach Program
tal_prog_file = fullfile(fileparts(which('gift.m')), 'icatb_talairach_scripts', 'talairach', 'talairach.jar');

talairach = fullfile(outputDir, [file_name, '_', flagCheck, '_', TALAIRACH_FILE, '.txt']);
talairach_z = fullfile(outputDir, [file_name, '_', flagCheck, '_', TALAIRACH_FILE, '_z.txt']);
talairach_excel = fullfile(outputDir, [file_name, '_', flagCheck, '_', TALAIRACH_EXCEL, '.xls']);
excel_file = [talairach_z, '.td'];

% Clean up files
cleanup_files({talairach, talairach_z, talairach_excel, excel_file});

check_ind = find(data > threshold);

if isempty(check_ind)
    disp(['No voxels found for this threshold (', num2str(threshold), ')']);
    return;
end

%The actual commands that will create the output files. Keep in mind
%that the second write_talairach command can give an error when there are no negative
%values to be found.
msgStr = 'Writing initial talairach files...';
disp(msgStr);
icatb_write_talairach(data, sample_size, threshold, xdim, ydim, zdim, voxel1, voxel2, voxel3, origin, talairach, talairach_z);

drawnow;
%This part will create the system command for execution.
msgStr = 'Running the Talairach Daemon... This part might take some time so please be patient.';
disp(msgStr);

commandStr = ['java -classpath "', tal_prog_file, '" org.talairach.ExcelToTD host=localhost:1600 2, "', talairach_z, '"'];

[status, result] = system(commandStr);

drawnow;
if (status == 1)
    error(result);
end

%Creation of the actual excel files.
C = readTextFile(excel_file, 8, '%n\t%n\t%n\t%s\t%s\t%s\t%s\t%s', '\t');

T = readTextFile(talairach_z, 5, '%n\t%n\t%n\t%n\t%s', '\t');

%If all goes well this command should make perfect sense.
C(1, 9) = T(1, 4);

C{1, 10} = T{1, 5};

%The script writes its ouput files in the directory you were working in.
disp('Creating Final Talairach Outputs in Excel Format');
icatb_generate_tal(C, talairach_excel, voxsize, volStr, valStr);
disp(['File is saved with the name ', talairach_excel]);
fprintf('\n');

% Clean up files
cleanup_files({talairach, excel_file});

function C = readTextFile(txtFile, numCols, formatStr, delimiter)
% Read text file using the format string

C = cell(1, numCols);

outStr = '';
for n = 1:numCols
    outStr = [outStr, 'C{', num2str(n), '}', ' '];
end

eval(['[', outStr, '] = textread(txtFile, formatStr, -1, ''delimiter'', delimiter);'])


function cleanup_files(files)
% Cleanup files

for nF = 1:length(files)
    current_file = files{nF};
    if exist(current_file, 'file')
        delete(current_file);
    end
end


function setJVM
% Set default JVM to be the one in matlabroot/sys/java

dirName = fullfile(matlabroot, 'sys', 'java', 'jre');

dirName = icatb_listSubFolders(dirName);

dirName = dirName(1, :);

d = dir(fullfile(dirName, 'jr*'));
d = d([d.isdir] == 1);

if ~isempty(d)
    dirName = fullfile(dirName, d(1).name);
    
    dirName = fullfile(dirName, 'bin');
    
    p = getenv('PATH');
    if ispc
        chk = findstr(lower(p), lower(dirName));
    else
        chk = findstr(p, dirName);
    end
    
    if (isempty(chk))
        setenv('PATH', [dirName, pathsep, p]);
    end
    
end


function [threshold, sample_size, convertToZ] = getParams(threshold, sample_size)
% Get inputs



% dialog Title
dlg_title = 'Talairach parameters';

numParameters = 1;

inputText(numParameters).promptString = 'Do you want to convert image to Z-scores?';
inputText(numParameters).uiType = 'popup';
inputText(numParameters).answerString = str2mat('No', 'Yes');
inputText(numParameters).dataType = 'string';
inputText(numParameters).tag = 'zscores';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;

% define all the input parameters in a structure
inputText(numParameters).promptString = 'Enter threshold';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(threshold);
inputText(numParameters).tag = 'threshold';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;

numParameters = numParameters + 1;

inputText(numParameters).promptString = 'Distance between contiguous voxels in mm';
inputText(numParameters).uiType = 'edit';
inputText(numParameters).answerString = num2str(sample_size);
inputText(numParameters).tag = 'dist';
inputText(numParameters).enable = 'on';
inputText(numParameters).value = 1;


% Input dialog box (get the necessary numbers)
answer = icatb_inputDialog('inputtext', inputText, 'Title', dlg_title, 'windowStyle', 'modal');


if ~isempty(answer)
    convertToZ = strcmpi(answer{1}, 'yes');
    threshold = str2num(answer{2});
    sample_size = str2num(answer{3});
else
    error('Figure window was quit');
end


function inputFiles = formFullPaths(inputFiles)
%% Form full paths

oldDir = pwd;

for nFile = 1:length(inputFiles)
    
    cF = inputFiles{nFile};
    [p, fN, extn] = fileparts(deblank(cF));
    
    if (isempty(p))
        p = pwd;
    end
    
    cd(p);
    
    inputFiles{nFile} = fullfile(pwd, [fN, extn]);
    
    cd(oldDir);
    
end