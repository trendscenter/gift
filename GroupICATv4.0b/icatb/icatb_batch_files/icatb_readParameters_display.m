function [parameters] = icatb_readParameters_display(inputFile)
% Function to read the parameters that will be used for displaying
% components.
%
% Inputs:
% inputFile - file containing display parameters. This must be a full file
% path.
%
% Output:
% parameters - Structure containing the display parameters

icatb_defaults;
global ANATOMICAL_PLANE;
global IMAGES_PER_FIGURE;
global COMPONENT_NAMING;

disp(['Reading parameters for displaying components ...']);

anatomicalPlane = ANATOMICAL_PLANE;
imagesPerFigure = str2num(IMAGES_PER_FIGURE);

parameters = struct;

% read display type
keywd = 'displayType';
inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'string');
% get the output directory field
displayType = getfield(inputData, keywd);
clear inputData;

if ~strcmpi(displayType, 'component explorer') & ~strcmpi(displayType, 'composite viewer') & ...
        ~strcmpi(displayType, 'orthogonal viewer')

    error('ErrorCheck:DisplayType', 'Not a valid option for displayType variable. \n See inputFile: %s for the options', ...
        inputFile);

end

parameters.displayType = displayType;

%%%%% Read input and output data %%%%%

%% Read data for ortho viewer
fmriFiles = [];

if strcmpi(displayType, 'orthogonal viewer')
    % Read input data
    keywd = 'sourceDir';
    inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'directory');
    sourceDir = getfield(inputData, keywd);
    clear inputData;

    % Read file pattern
    keywd = 'sourceFilePattern';
    inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'string');
    sourceFilePattern = getfield(inputData, keywd);
    clear inputData;

    fmriFiles = icatb_listFiles_inDir(sourceDir, sourceFilePattern);

    if isempty(fmriFiles)
        error(['Please check the file pattern for fmri data in variable sourceFilePattern']);
    end

    fmriFiles = icatb_fullFile('directory', sourceDir, 'files', fmriFiles);

end

parameters.fmriFiles = fmriFiles;


%% Read component data
keywd = 'outputDir';
inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'directory');
outputDir = getfield(inputData, keywd);
clear inputData;

parameters.outputDir = outputDir;
parameters.filesOutputDir = parameters.outputDir;

% Read file pattern
keywd = 'compFilePattern';
inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'string');
compFilePattern = getfield(inputData, keywd);
clear inputData;

compFiles = icatb_listFiles_inDir(outputDir, compFilePattern);
if isempty(compFiles)
    error(['Please check the file pattern for components data in variable compFilePattern']);
end

checkCompFile = find(icatb_regexpm(cellstr(compFiles), COMPONENT_NAMING));

if length(checkCompFile) ~= size(compFiles, 1)
    error(['Component file pattern must contain ', COMPONENT_NAMING, ' pattern']);
end

compFiles = icatb_fullFile('directory', outputDir, 'files', compFiles);

parameters.compFiles = compFiles;

compNumbers = [];
keywd = 'compNumbers';
try
    inputData = icatb_read_variables(inputFile, keywd, 'vector', 'integer');
    compNumbers = getfield(inputData, keywd);
catch
end

if min(compNumbers) == 0
    error('compNumbers variable cannot have value of zero');
end

clear inputData;


% Count number of components
[allComp] = icatb_get_countTimePoints(compFiles);

if isempty(compNumbers)
    compNumbers = [1:1:allComp];
end

if max(compNumbers) > allComp
    error(['compNumbers variable maximum exceeds the number of components']);
end

if strcmpi(displayType, 'orthogonal viewer')
    if length(compNumbers ) > 1
        disp(['Can display only one component using orthogonal viewer ...']);
        compNumbers = compNumbers(1);
    end
end

if strcmpi(displayType, 'composite viewer')
    if length(compNumbers) > 5
        disp(['Can display only atmost five components using composite viewer ...']);
        compNumbers = compNumbers(1:5);
    end
end

disp(['Selected component/components is/are: ', num2str(compNumbers)]);

parameters.compNumbers = compNumbers;
parameters.numComp = length(parameters.compNumbers);

keywd = 'structFile';
inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'file');
structFile = getfield(inputData, keywd);
clear inputData;

parameters.structFile = structFile;

%%% Done reading input and output data %%%%%


%%%% Display parameters %%%%%%%%

keywd = 'returnValue';
inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'integer');
returnValue = getfield(inputData, keywd);
clear inputData;

if returnValue > 4
    returnValue = 4;
end

parameters.returnValue = returnValue;
parameters.imagevalues = returnValue;


keywd = 'convertToZ';
inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'integer');
convertToZ = getfield(inputData, keywd);
clear inputData;

if convertToZ ~= 0
    convertToZ = 1;
end

parameters.convertToZ = convertToZ;

keywd = 'thresholdValue';
inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'numeric');
thresholdValue = getfield(inputData, keywd);
clear inputData;

thresholdValue = abs(thresholdValue);

parameters.thresholdvalue = thresholdValue;

if strcmpi(displayType, 'component explorer')
    keywd = 'imagesPerFigure';
    inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'integer');
    imagesPerFigure = getfield(inputData, keywd);
    clear inputData;
end

parameters.imagesperfigure = imagesPerFigure;

if ~strcmpi(displayType, 'orthogonal viewer')
    keywd = 'anatomicalPlane';
    inputData = icatb_read_variables(inputFile, keywd, 'scalar', 'string');
    anatomicalPlane = getfield(inputData, keywd);
    clear inputData;

    if ~strcmpi(anatomicalPlane, 'axial') & ~strcmpi(anatomicalPlane, 'sagittal') & ...
            ~strcmpi(anatomicalPlane, 'coronal')
        error('anatomicalPlane variable must be axial or saggital or coronal');
    end

end

parameters.anatomicalplane = anatomicalPlane;

%%%%%% End for display parameters %%%

parameters.imagevalues = parameters.returnValue;
tmp = icatb_get_vol_nifti(parameters.structFile);
tmp = icatb_get_slice_def(tmp, parameters.anatomicalplane);
parameters.slicerange = tmp.slices;


disp(['Done reading']);