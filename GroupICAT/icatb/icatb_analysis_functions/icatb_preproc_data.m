function out = icatb_preproc_data(out, preProcType, verbose)
%% Pre-process data
%
% Inputs:
% 1. out - 2D Array (Voxels X time)
% 2. preProcType:
%    1 - Remove mean per time point
%    2 - Remove mean per voxel
%    3 - Intensity normalization
%    4 - Variance normalization
%

options = {'Remove Mean Per Timepoint', 'Remove Mean Per Voxel', 'Intensity Normalization', 'Variance Normalization'};
alias_options = {'rt', 'rv', 'in', 'vn'};
if (nargin == 0)
    out = options;
    return;
end

if (~exist('verbose', 'var'))
    verbose = 1;
end

if (isnumeric(preProcType))
    preProcType = options{preProcType};
end

preProcType = lower(preProcType);

ind = strmatch(preProcType, lower(options), 'exact');
if (isempty(ind))
    ind = strmatch(preProcType, lower(alias_options), 'exact');
end
msgStr = options{ind(1)};

%% Use 50 voxels at a time for variance or intensity normalization
blockSize = 50;
voxels = size(out, 1);
nLoops = ceil(voxels / blockSize);

msg = '';
switch (preProcType)
    case {'remove mean per timepoint', 'rt'}
        msg = 'Removing mean per time point ...';
        %% Remove mean per time point
        out = icatb_remove_mean(out);
    case {'remove mean per voxel', 'rv'}
        msg = 'Removing mean per voxel ...';
        %% Remove mean per voxel
        out = icatb_remove_mean(out')';
    case {'intensity normalization', 'variance normalization', 'in', 'vn'}
        %% Intensity or Variance normalization
        
        msg = ['Using ', msgStr, ' ...'];
        endT = 0;
        
        %% Use blocks
        for n = 1:nLoops
            
            startT = endT + 1;
            endT = endT + blockSize;
            
            if (endT > voxels)
                endT = voxels;
            end
            
            tmp = out(startT:endT, :)';
            
            if (strcmpi(preProcType, 'intensity normalization') || strcmpi(preProcType, 'in'))
                %% Intensity normalization
                tmp = repmat(100./(mean(tmp) + eps), size(tmp, 1), 1) .* tmp;
                
                out(startT:endT, :) = tmp';
                
                clear tmp;
                
            else
                %% Variance normalization
                
                tmp = detrend(tmp);
                tmp = tmp.*repmat(1./(std(tmp) + eps), size(tmp, 1), 1);
                
                out(startT:endT, :) = tmp';
                
                clear tmp;
                
            end
            
        end
        
end

if (verbose)
    disp(msg);
end