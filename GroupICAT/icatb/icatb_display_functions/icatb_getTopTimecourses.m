function  [A, values, topComponents] = icatb_getTopTimecourses(Voxel, ICAFiles, structFile, dataType, complexInfo, ...
zipFileName, files_in_zip)
% gets the top five contibutors for a particular voxel
%
% Input:
% 1. Voxel - Voxel location
% 2. ICAFiles - ICA files
% 3. StructFile - structural file
%
% Output:
% 1. A - Time course
% 2. values - voxel values
% 3. topComponents - top five contributors

if ~exist('dataType', 'var')
    dataType = 'real';
end

if ~exist('complexInfo', 'var')
    complexInfo = [];
end

if ~exist('zipFileName', 'var')
    zipFileName = {};
    files_in_zip = {};
end

icatb_defaults;
global DETRENDNUMBER;

numberOfTopComp = 5;

% ICA component files
P = ICAFiles;

% get the voxel values
[voxelValues] = icatb_getVoxelValues(Voxel, P, structFile, dataType, complexInfo, 'write', zipFileName, files_in_zip);

% component indices
%componentIndex = [1:size(P,1)];
componentIndex = [1:length(voxelValues)];

% get the sorted values
[values indices] = sort(abs(voxelValues));

% top component indices
indices = indices(end:-1:1);
topComponents = indices(1:numberOfTopComp);
values = voxelValues(topComponents);

% get the ICA time course
A = icatb_loadICATimeCourse(ICAFiles, dataType, complexInfo, topComponents);

% get only the real or magnitude part
if isa(A, 'complex_data')
    A = getfield(A, 'firstField');
end

% % detrend ICA time course
% A = icatb_detrend(A, 1, size(A, 1), DETRENDNUMBER);
