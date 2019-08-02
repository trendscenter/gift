function icatb_save(varargin)
%% Function to save the MAT files. The function inputs are the same as save command. For versions greater than MATLAB 6.5,
% option is provided to save the MAT files in the appropriate version.
%

%% Load defaults
icatb_defaults;

%% Enforce MAT file version
global ENFORCE_MAT_FILE_VER;

%% MATLAB version
matlab_version = icatb_get_matlab_version;

%% Loop over nargs
for ii = 2:nargin
    if ~strcmp(varargin{ii}(1), '-')
        % get only the variable names
        eval([varargin{ii}, ' = ', 'evalin(''caller'', varargin{ii});']);
    end
end

p = fileparts(varargin{1});
if (~isempty(p))
    if (exist(p, 'dir') ~= 7)
        [parDir, childDir] = fileparts(p);
        mkdir(parDir, childDir);
    end
end

%% Enforce MAT file version
if (matlab_version <= 13)
    ENFORCE_MAT_FILE_VER = '';
end

if (~isempty(ENFORCE_MAT_FILE_VER))
    save(varargin{1}, varargin{2:end}, ENFORCE_MAT_FILE_VER);
else
    save(varargin{1}, varargin{2:end});
end