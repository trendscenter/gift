function [A] = icatb_loadICATimeCourse(P, dataType, complexInfo, compIndicies, ...
    zipFileName, files_in_zip)
% Loads ICA Time course
%
% Input:
% 1. P - component images naming
% 2. dataType - 'real' or 'complex'
% 3. complexInfo - structure containing image namings and complex type
% ('real&imaginary' or 'magnitude&phase')
% Output:
% A - time course of class type double or complex_data classd

if ~exist('P', 'var')
    P = icatb_selectEntry('typeSelection', 'multiple', 'filter', '*comp*.img', 'title', 'Select Component file/files', 'filetype', 'image', ...
        'filenumbers', 1);
end

if ~exist('dataType', 'var')
    dataType = 'real';
end

if ~exist('zipFileName', 'var')
    zipFileName = [];
end

if ~exist('files_in_zip', 'var')
    files_in_zip = {};
end

if ~exist('compIndicies', 'var')
    compIndicies = [];
end

icatb_defaults;

global COMPONENT_NAMING;
global TIMECOURSE_NAMING;
global FUNCTIONAL_DATA_FILTER;
% variables to detect the reading and writing of complex images
global READ_NAMING_COMPLEX_IMAGES;
global WRITE_NAMING_COMPLEX_IMAGES;

modalityType = icatb_get_modality;

P = icatb_parseExtn(deblank(P(1, :)));

% image extension
[pathstr_comp, file_comp, imExtn] = fileparts(P);

if ~isempty(zipFileName)
    icatb_unzip(regexprep(zipFileName, ['.*\', filesep], ''), pathstr_comp);
end

lastUnderScore = icatb_findstr(deblank(P(1,:)),'_');
lastUnderScore = lastUnderScore(end);
component_name = deblank(P(1,1:lastUnderScore));
timecourse_name = strrep(component_name, COMPONENT_NAMING, TIMECOURSE_NAMING);

file = [timecourse_name, imExtn];

% get the file parts
[pathstr, fileName, extn] = fileparts(file);

% get the present working directory
if isempty(pathstr)
    pathstr = pwd;
end

% form full file name
fileName = fullfile(pathstr, [fileName, imExtn]);


% check the data type and load timecourses accordingly
if strcmpi(dataType, 'complex')

    if ~exist('complexInfo', 'var')
        error('Complex Info is not present');
    end

    %     % complex file naming
    complex_file_naming = complexInfo.complex_file_naming;
    % complex type
    complex_type = complexInfo.complexType;

    % get the file namings
    [P] = icatb_get_complex_files_naming(fileName, dataType, complexInfo, 'write');

    % form naming for first set
    PFirst = str2mat(P.first);

    % read data and header info of first set
    [data] = icatb_loadData(PFirst);

    % form naming for second set
    PSecond = str2mat(P.second);

    % read data and header info of second set
    [data2] = icatb_loadData(PSecond);

    % get the component indicies
    if isempty(compIndicies)
        compIndicies = (1:1:size(data2, 2));
    end

    %     % get the corresponding time course
    firstField = zeros(size(data, 1), length(compIndicies));
    secondField = firstField;
    % loop over selected components
    for ii = 1:length(compIndicies)
        firstField(:, ii) = data(:, compIndicies(ii));
        secondField(:, ii) = data2(:, compIndicies(ii));
    end

    clear data; clear data2;

    % form time course structure
    A = struct('firstField', firstField, 'secondField', secondField);

    % form object of class complex_data
    A = complex_data(A);

else
    % load real images

    % get the volume information
    allTC = icatb_loadData(fileName);

    % get the component indicies
    if isempty(compIndicies)
        compIndicies = (1:1:size(allTC, 2));
    end

    %     % Initialise time course
    A = zeros(size(allTC, 1), length(compIndicies));
    % get the corresponding time course
    for ii = 1:length(compIndicies)
        A(:, ii) = allTC(:, compIndicies(ii));
    end

end
% end for checking data type

% delete the files unzipped
if ~isempty(zipFileName)
    files_in_zip = str2mat(regexprep(files_in_zip, ['.*\', filesep], ''));
    %fileN = icatb_fullFile('files', files_in_zip, 'directory', pathstr_comp);
    icatb_delete_file_pattern(files_in_zip, pathstr_comp);
end