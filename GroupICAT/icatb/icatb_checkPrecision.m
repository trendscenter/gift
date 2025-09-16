function precisionType = icatb_checkPrecision(precisionType, verbose)
%% Check precision
%

if (~exist('verbose', 'var'))
    verbose = 1;
end

precisionType = lower(precisionType);

% Matlab version
matlab_version = icatb_get_matlab_version;

% Get mex extension
mexExtension = mexext;

% Operating System bit
OSBIT = 32;
if (strcmp(mexExtension(end-1:end), '64'))
    OSBIT = 64;
end

%% Convert matrix to double taking MATLAB version and Operating System bit
% into account
if ((matlab_version == 13) && (strcmpi(precisionType, 'single')))
    if (verbose)
        disp('Single precision is not supported on R13. Using double instead ...');
    end
    precisionType = 'double';
end

if ((OSBIT == 64) && (matlab_version < 2008) && (strcmpi(precisionType, 'single')))
    if (verbose)
        disp('Single precision won''t work on MATLAB versions < R2008a for 64 bit OS. Using double instead ...');
    end
    precisionType = 'double';
end