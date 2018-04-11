function icatb_componentExplore(parameters, minimumParametersStructure)

% component explorer displays components of a particular viewing set.
% This method can be used by using display GUI or by selecting
% component explorer button.
% Input: parameters containing the following fields
% 1. returnValue - return value of the image values ('Positive and Negative' means 1, 'Positive',
% means 2,  'Absolute Value' means 3 and 'Negative' means 4).
% 2. anatomicalplane - 'Axial', 'Sagital' or 'Coronal'
% 3. imagesperfigure - '1', '4', '16', etc.
% 4. slicerange - minimum:interval:maximum
% 5. numComp - number of components to be displayed
% 6. thresholdvalue - given threshold
% 7. structuralImage - 3D matrix
% 8. structHInfo - Header information
% 9. convertToZ - converttozscores value ('Yes' means 1, 'No' means 0)
% 10. icasig - component image as 3D matrix
% 11. sortParameters - structure contains the sorting information (see
% icatb_sortComponentsGUI and icatb_sortComponents)
%
% Output: displays components using icatb_drawMComponents

% Load display parameters window to get input from the user
if ~exist('parameters', 'var') | exist( 'minimumParametersStructure', 'var' )

    icatb_defaults;
    global DETRENDNUMBER;
    global SMOOTHPARA;
    global SMOOTHINGVALUE;

    files_in_zip = {};

    if (exist( 'minimumParametersStructure', 'var' ) & ~isempty(  minimumParametersStructure ))
        parameters = minimumParametersStructure;
        clear minimumParametersStructure;
    else
        % load the component explorer GUI
        parameters = icatb_inputParameters_display('component explorer');
        if isempty(parameters)
            error('figure window was quit');
        end
    end

    imageValues = parameters.returnValue; % image values
    slicePlane = parameters.anatomicalplane; % slice plane
    sliceRange = parameters.slicerange; % slice range in mm
    numComp = parameters.numComp; % Number of components
    threshValue = parameters.thresholdvalue; % threshold value
    compFiles = parameters.compFiles; % component files
    compNumbers = parameters.compNumbers; % component numbers
    convertToZ = parameters.convertToZ; % convert to z scores
    structuralFile = parameters.structFile; % structural file

    if isfield(parameters, 'files_in_zip')
        % get files in zip
        files_in_zip = parameters.files_in_zip;
    end

    % load the component files here
    [icasig, A, structuralImage, coords, HInfo, parameters.text_left_right] = icatb_loadICAData('structFile', ...
        structuralFile, 'compFiles', compFiles, 'slicePlane', slicePlane, 'sliceRange', sliceRange, 'comp_numbers', ...
        compNumbers, 'convertToZ', convertToZ, 'returnValue', imageValues, 'threshValue', ...
        threshValue, 'dataType', 'real', 'complexInfo', []);
    %%%%%%%%% Apply structural image parameters %%%%%%%%


    % delete the unzipped files
    if ~isempty(files_in_zip)
        for ii = 1:length(files_in_zip)
            drawnow;
            delete(files_in_zip{ii});
        end
    end


    % smooth ica time course
    if strcmpi(SMOOTHPARA, 'yes')
        icaTimecourse = icatb_gauss_smooth1D(icaTimecourse, SMOOTHINGVALUE);
    end

    structHInfo = HInfo;
    parameters.structHInfo = HInfo;
    structDIM = HInfo.DIM;
    parameters.structuralImage = structuralImage;
    parameters.undetrendICA = A; % show the user detrend and undetrend
    parameters.icaTimecourse = icatb_detrend(A, 1, size(A, 1), DETRENDNUMBER); % detrended time course
    parameters.modelTimecourse = [];
    parameters.htmlFile = 'icatb_component_explorer.htm';
    [pp, viewingSetName, extn] = fileparts(deblank(compFiles(1, :)));
    parameters.figLabel = viewingSetName;
    for ii = 1:numComp
        compLabels(ii).string = ['Component ', num2str(compNumbers(ii))];
    end
    parameters.compLabels = compLabels;
    clear prefix;
    % specify a flag here that the component explorer is accessed
    % independently
    flagComponentExplorer = 'no_displayGUI';

else
    imageValues = parameters.imagevalues;  % image values
    slicePlane = parameters.anatomicalplane; % slice plane
    sliceRange = parameters.slicerange; % slice range
    plotCount = parameters.imagesperfigure; % images per figure
    numComp = parameters.numComp;
    threshValue = parameters.thresholdvalue;
    structuralImage = parameters.structuralImage;
    parameters = rmfield(parameters, 'structuralImage'); % remove the field structural
    structHInfo = parameters.structHInfo;
    convertToZ = parameters.convertToZ;
    icasig = parameters.icasig;  % get the component images
    parameters = rmfield(parameters, 'icasig'); % remove images from parameters structure
    otherFields = {'modelTimecourse', 'refInfo', 'diffTimePoints'}; % fields to be added
    for ii = 1:length(otherFields)
        if ~isfield(parameters, otherFields{ii})
            parameters = setfield(parameters, otherFields{ii}, []);
        end
    end
    % specify a flag here that the component explorer is accessed
    % using display GUI
    flagComponentExplorer = 'displayGUI';
end

% defaults
icatb_defaults;
global COLORMAP_FILE;
load(COLORMAP_FILE);

% %get colormap
cm = icatb_getColormap(1, imageValues, 1);

% return overlayed images depending upon the data type of icasig and
% structuralImage
[images,  maxICAIM, minICAIM, minInterval, maxInterval] = icatb_check_overlayComp(icasig, structuralImage, ...
    imageValues, numComp);

size_structuralImage = size(structuralImage);
clear icasig; clear im;
clear structuralImage;

%%% Put these additional checks for single slice in sagittal and coronal planes
if strcmpi(slicePlane, 'sagittal')
    if length(size_structuralImage) == 2
        size_structuralImage = [1, size_structuralImage];
    end
end

if strcmpi(slicePlane, 'coronal')
    if length(size_structuralImage) == 2
        size_structuralImage = [size_structuralImage(1), 1, size_structuralImage(2)];
    end
end
%%%%%% End for putting slice plane checks %%%%%

% reshape images according to data type
reshapeSize = [numComp, size_structuralImage];
images = reshape(images, reshapeSize);


%--get image in correct plane
if strcmpi(slicePlane, 'sagittal')
    images = permute(images, [1 3 4 2]);
    %structHInfo.DIM = [structHInfo.DIM(2) structHInfo.DIM(3) structHInfo.DIM(1)];
    structHInfo.DIM = [size_structuralImage(2), size_structuralImage(3), size_structuralImage(1)];
    structHInfo.VOX = [structHInfo.VOX(2) structHInfo.VOX(3) structHInfo.VOX(1)];
end

if strcmpi(slicePlane, 'coronal')
    images = permute(images, [1 2 4 3]);
    %images = icatb_permuteImage(images, [1 2 4 3]);
    %structHInfo.DIM = [structHInfo.DIM(1) structHInfo.DIM(3) structHInfo.DIM(2)];
    structHInfo.DIM = [size_structuralImage(1), size_structuralImage(3), size_structuralImage(2)];
    structHInfo.VOX = [structHInfo.VOX(1) structHInfo.VOX(3) structHInfo.VOX(2)];
end

% check for the number of data sets concatenated
if ~isfield(parameters, 'numSubjects')
    parameters.numSubjects = 1;
end
if ~isfield(parameters, 'numSessions')
    parameters.numSessions = 1;
end
parameters.num_DataSets = parameters.numSubjects*parameters.numSessions;

% append to parameters structure
parameters.images = images; parameters.cm = cm;
parameters.minICAIM = minICAIM; parameters.maxICAIM = maxICAIM;
parameters.minInterval = minInterval; parameters.maxInterval = maxInterval;

% % cleaned drawMcomponents
icatb_drawMComponents(parameters);