function files = icatb_rename_4d_file(files, fileNumber)
%For nifti and analyze files rename the files by adding a number at the
%end

if ~exist('fileNumber', 'var')
    fileNumber = [];
end

if ~isempty(files)
    
    % Do pattern match to check IMG or NII files
    files = cellstr(files);
    
    % MATCH IMG OR NII
    checkNII = regexpi(files, '\.nii$|\.gz$');
    checkIMG = regexpi(files, '\.img$');
    
    good_inds1 = icatb_good_cells(checkNII);
    good_inds2 = icatb_good_cells(checkIMG);
    
    % Get good cells
    good_inds = (good_inds1 | good_inds2);
    good_inds = find(good_inds ~= 0);
    
    clear good_inds1 good_inds2 checkNII checkIMG;
    
    % If IMG or NII files exist check the headers of the image files
    if ~isempty(good_inds)
        imFiles = files(good_inds);
        if (~isempty(fileNumber)) && (length(fileNumber) == 1) && (fileNumber == 1)
            imFiles = strcat(imFiles, ',1');
            files(good_inds) = imFiles;
        else
            % Loop over img or nii files
            for nIm = 1:length(imFiles)
                currentFile = imFiles{nIm};
                
                try
                    numFiles = getDims(currentFile);
                catch
                    numFiles = 1;
                end
                
                if (isempty(fileNumber))
                    tempFileNum = (1:numFiles);
                else
                    tempFileNum = fileNumber;
                end
                
                tempFileNum(tempFileNum > numFiles) = [];
                
                if ~isempty(tempFileNum)
                    tempFiles = [repmat(currentFile, length(tempFileNum), 1), repmat(',', length(tempFileNum), 1), numberToString(tempFileNum(:))];
                else
                    tempFiles = '';
                end
                files{good_inds(nIm)} = tempFiles;
            end
            % End loop over img or nii files
            files = cellstr(char(files));
            ind = icatb_good_cells(files);
            files = files(ind);
        end
    end
    % End for checking IMG or NII files
    
    if ~isempty(files)
        files = char(files);
    end
    
end


function str = numberToString(nums)
% number to string. pad space after the number

try
    str = arrayfun(@num2str, nums, 'uniformoutput', false);
catch
    str = cell(length(nums), 1);
    for n = 1:length(nums)
        str{n} = num2str(nums(n));
    end
end

str = char(str);


function tp = getDims(currentFile)

if (~strcmpi(currentFile(end-2:end),'.gz'))
ni = icatb_nifti(currentFile);
tp = ni.dat.dim(4);
else
    ni = icatb_read_hdr(currentFile);
    tp = ni.dim(5);
end