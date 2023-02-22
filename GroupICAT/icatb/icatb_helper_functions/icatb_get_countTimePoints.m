function [diffTimePoints, extns, dims] = icatb_get_countTimePoints(files)
% Returns the number of time points in files
%
% Input:
% files - array of structures or character array in columns (output of
% str2mat function)
%
% Output:
%

dims = [];
if ~isstruct(files)
    % check if the files are in character format
    if ischar(files)
        temp = files;
        clear files;
        files.name = temp;
        clear temp;
    else
        error('Files should be in structure or character format');
    end
end

checkDim = 0;
diffTimePoints = zeros(1, length(files));
% loop over files
for ii = 1:length(files)
    currentFiles = files(ii).name; % current subject file
    countTimePoints = 0;
    % loop over the number of image files
    for jj = 1:size(currentFiles, 1)
        checkDim = checkDim + 1;
        [pathstr, fName, extn] = fileparts(deblank(currentFiles(jj, :)));
        if isempty(pathstr)
            pathstr = pwd;
        end
        % store the extension
        extns.data(ii).file(jj).ext = extn;
        extns.data(ii).file(jj).file_name = fName;
        extns.data(ii).file(jj).dir = pathstr;
        % get the image extension
        if strcmpi(extns.data(ii).file(jj).ext, '.nii') || strcmpi(extns.data(ii).file(jj).ext, '.gz')
            % read the time point dimensions
            [hdr] = icatb_read_hdr(deblank(currentFiles(jj, :)));
            % update the time points information
            try
                countTimePoints = countTimePoints + hdr.dime.dim(5);
                dims = hdr.dime.dim(2:4);
            catch
                countTimePoints = countTimePoints + hdr.dim(5);
                dims = hdr.dim(2:4);
            end            
        elseif strcmpi(extns.data(ii).file(jj).ext, '.img')
            % update count time points
            %[pathstr, hdrFile, extn] = fileparts(deblank(currentFiles(jj, :)));
            hdrFile = fName;
            hdrFile = fullfile(pathstr, [hdrFile, '.hdr']);
            if ~exist(hdrFile, 'file')
                error(['Header file ', hdrFile, ' doesn''t exist']);
            end
            [hdr, otherendian] = icatb_read_hdr(hdrFile);
            countTp = hdr.dime.dim(5);
            countTimePoints = countTimePoints + countTp;
            dims = hdr.dime.dim(2:4);
        else
            countTimePoints = countTimePoints + 1;
        end
        
        
        % get the dimension from the first file
        if checkDim == 1
            oldDim = dims;
        end
        
        % checking the dimensions
        if ~isempty(dims)
            
            checkVec = find((oldDim == dims) > 0);
            
            % check the data dimensions
            if length(checkVec) ~= length(oldDim)
                error(['Data dimensions should be the same for all files. Please check the dimension for file ', ...
                    deblank(currentFiles(jj, :))]);
            end
            % end for checking the data dimensions
            
        end
        % end for checking the dimensions
    end
    % end loop over the number of image files
    % store the time points information
    diffTimePoints(ii) = countTimePoints;
    % end loop over the files of the current subject
    
end
% end loop over all subject files
