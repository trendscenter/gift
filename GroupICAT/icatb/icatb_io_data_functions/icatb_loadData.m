function [data, HInfo, XYZ] = icatb_loadData(P, dataType, complexInfo, optional, file_number, mask_ind)
%% Function to load the data. Presently, 3D analyze, 4D analyze and Nifti
% data are handled.
% Uses spm_vol and spm_read_vols to get the data.
%
% Input:
%
% 1. P - full file path for the variables
% If complex data is to be read then other variables are needed
% 2. dataType - 'real' or 'complex'
% 3. complexInfo - Variable containing complex images prefix and complex
% type ('real&imaginary' or 'magnitude&phase')
% 4. optional - By default complex images are assumed to be of 'read' type.
% This variable allows only prefix for 'write' but allows both prefix and
% suffix for read type.
% 5. file_number - By default all images are loaded. file_number must be a
% vector of numbers.
%
% Ouput:
%
% 1. data - Real or complex depending upon the dataType variable
% 2. HInfo - Structure containing the volume information.

% by default data type is assumed to be real
if ~exist('dataType', 'var')
    dataType = 'real';
end

% by default assuming complex file naming to be of read type
if ~exist('optional', 'var')
    optional = 'read';
end

% file number
if ~exist('file_number', 'var')
    file_number = [];
end

if ~exist('mask_ind', 'var')
    mask_ind = [];
end

% check for empty args
if isempty(dataType)
    dataType = 'real';
end

%% Check real or imaginary
if strcmpi(dataType, 'real')

    % data and header dimensions
    [data, HInfo, XYZ] = icatb_read_data(P, file_number, mask_ind);

else
    % if data type is complex

    if ~exist('complexInfo', 'var')
        error('Complex Info variable is not passed.');
    end

    % complex type
    complex_type = complexInfo.complexType;

    [P] = icatb_get_complex_files_naming(P, dataType, complexInfo, optional);

    % form naming for first set
    PFirst = str2mat(P.first);

    % read data and header info of first set
    [data, HInfo, XYZ] = icatb_read_data(PFirst, file_number, mask_ind);

    % form naming for second set
    PSecond = str2mat(P.second);

    % read data and header info of second set
    [data2, HInfo2, XYZ] = icatb_read_data(PSecond, file_number, mask_ind);

    % depending upon the the complex type form complex data
    if strcmpi(complex_type, 'real&imaginary')
        % deal real and imaginary case
        data = complex(data, data2);
        clear data2;
    else
        % deal magnitude and phase case
        data = complex(data.*cos(data2), data.*sin(data2));
        clear data2;
    end
    % end for complex type

end
% end for checking real or imaginary
