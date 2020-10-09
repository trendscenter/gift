function icatb_parSave(fname, vars, varnames, optional)
%% Save variables in parallel mode
%

%% Load defaults
icatb_defaults;

% Enforce MAT file version
global ENFORCE_MAT_FILE_VER;


if (~exist('optional', 'var'))
    optional = '';
end

% MATLAB version
matlab_version = icatb_get_matlab_version;


for i = 1:length(vars)
    eval([varnames{i}, '=vars{i};']);
end

p = fileparts(fname);

if (~isempty(p))
    if (exist(p, 'dir') ~= 7)
        [parDir, childDir] = fileparts(p);
        mkdir(parDir, childDir);
    end
end

% Enforce MAT file version
if (matlab_version <= 13)
    ENFORCE_MAT_FILE_VER = '';
end

enforce_mat_ver = ENFORCE_MAT_FILE_VER;

%% Append flags if any
if (~isempty(optional))
    varnames{end + 1} = optional;
    if (strcmpi(optional, '-append'))
        enforce_mat_ver = '';
    end
end

if (~isempty(enforce_mat_ver))
    varnames{end + 1} = enforce_mat_ver;
end

%% Save variables
save(fname, varnames{:});