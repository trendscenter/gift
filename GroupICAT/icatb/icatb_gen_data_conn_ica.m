function sesInfo = icatb_gen_data_conn_ica(sesInfo)
%% Generate connectivity matrices and save the files in MAT format
%

outputDir = sesInfo.userInput.pwd;
prefix = sesInfo.userInput.prefix;
filesInfo = sesInfo.userInput.dataInfo.filesInfo;
numOfSess = filesInfo.numOfSess;
allFiles = cellstr(filesInfo.filesList);
numOfDataSets = length(allFiles);
numOfSub = ceil(numOfDataSets / numOfSess);
dummy_scans = filesInfo.file_numbers;
if (isempty(dummy_scans))
    dummy_scans = 0;
end

conn_dir = 'connectivity_matrices';

maskFile = sesInfo.userInput.dataInfo.maskFile;

%% Check the filepatterns and generate file structure

files = repmat(struct('name', []), 1, numOfSub*numOfSess);
diffTimePoints = zeros(1, numOfSub*numOfSess);

input_data_file_patterns = cell(numOfSub, numOfSess);
endTp = 0;
for nSub = 1:numOfSub
    startTp = endTp + 1;
    endTp = endTp + numOfSess;
    input_data_file_patterns(nSub, :) = allFiles(startTp:endTp);
end

counter = 0;
for nSub = 1:size(input_data_file_patterns, 1)
    for nSess = 1:size(input_data_file_patterns, 2)
        counter = counter + 1;
        current_file_name = input_data_file_patterns{nSub, nSess};
        [pathstr, fileN, extn] = fileparts(current_file_name);
        fileList = icatb_listFiles_inDir(pathstr, [fileN, extn]);
        fileListWithDir = icatb_fullFile('files', fileList, 'directory', pathstr);
        
        fileListWithDir = icatb_rename_4d_file(fileListWithDir);
        
        if (length(dummy_scans) > 1)
            fileNum = dummy_scans;
            fileNum(fileNum > size(fileListWithDir, 1)) = [];
            if isempty(fileNum)
                error('Error:DummyScans', 'Cannot find the files specified with file numbers. Please check dummy_scans variable for \nfile %s\n', current_file_name);
            end
            fileListWithDir = fileListWithDir(fileNum, :);
        else
            if (dummy_scans > 0)
                if (dummy_scans >= size(fileListWithDir, 1))
                    error('Error:DummyScans', 'Please check dummy_scans variable (%d) as it exceeds or equals the no. of time points (%d) for\nfile %s\n',  ...
                        dummy_scans, size(fileListWithDir, 1), current_file_name);
                end
                fileNum = (1:size(fileListWithDir, 1));
                fileNum(1:dummy_scans) = [];
                fileListWithDir = fileListWithDir(fileNum, :);
                
            end
        end
        
        files(counter).name = fileListWithDir; % append the file list with directory
        diffTimePoints(counter) = size(fileListWithDir, 1);
    end
end

chk_sessInfo = sesInfo;
chk_sessInfo.userInput.dataType = 'real';
chk_sessInfo.userInput.files = files;
chk_sessInfo.userInput.maskFile = maskFile;
chk_sessInfo = icatb_update_mask(chk_sessInfo);

tmp_mask_name = fullfile(outputDir, [prefix, 'Mask.nii']);
maskV = icatb_spm_vol(tmp_mask_name);
maskV = maskV(1);
mask = icatb_spm_read_vols(maskV);
mask(isfinite(mask) == 0) = 0;

sesInfo.userInput.HInfo = chk_sessInfo.userInput.HInfo;

%% Subsample mask
subsampling_depth = sesInfo.userInput.dataInfo.subsampling_depth;

if (subsampling_depth > 1)
    
    new_mask = zeros(size(mask));
    new_mask(1:subsampling_depth:end, 1:subsampling_depth:end, 1:subsampling_depth:end) = 1;
    
    mask2 = new_mask.*mask;
    mask2 = double(mask2 ~= 0);
    
    
else
    
    mask2 = mask;
    
end

mask_ind = find(mask2 ~= 0);

num_voxels = length(mask_ind);
if (num_voxels < 10)
    error('Number of voxels found is less than 10. Check the subsampling depth or the mask');
end

tmp_mask_name = fullfile(outputDir, [prefix, 'Mask_subsampled.nii']);
maskV.fname = tmp_mask_name;
maskV.n(1) = 1;
icatb_write_vol(maskV, mask2);

sesInfo.userInput.mask_ind  = mask_ind;
sesInfo.userInput.inputFiles = files;
sesInfo.userInput.maskFile = tmp_mask_name;




%% Cleanup temporary batch file
tmp_batch = fullfile(outputDir, [prefix,  '_gift_batch_file.m']);
try
    if exist(tmp_batch, 'file')
        delete(tmp_batch);
    end
catch
end


%% Generate functional connectivity matrices
conn_type = 'ENLwFC';
try
    conn_type = sesInfo.userInput.dataInfo.conn_type;
catch
end

conn_input = fullfile(outputDir, conn_dir);

if (~exist(conn_input, 'dir'))
    mkdir(outputDir, conn_dir);
end

clear files


%sesInfo.userInput.HInfo = ;
mask_ind = find(mask2 ~= 0);

files = repmat(struct('name', []), 1, numOfSub*numOfSess);
parfor nDataset = 1:numOfDataSets
    
    subNum = ceil(nDataset/numOfSess);
    sesNum = mod(nDataset-1, numOfSess) + 1;
    
    
    disp(['Computing connectivity matrices of subject ', num2str(subNum), ' session ', num2str(sesNum), ' ...']);
    
    tmpFiles = sesInfo.userInput.inputFiles(nDataset).name;
    tmpDat = icatb_read_data(tmpFiles, [], mask_ind);
    tmpDat = tmpDat';
    conn_matrix = icatb_calc_ENLwFC(single(tmpDat));

    outFile = fullfile(outputDir, conn_dir, [prefix, '_sub', icatb_returnFileIndex(subNum), '_s', num2str(sesNum), '_FULL_conn_data.mat']);
    disp(['Saving file ', outFile, ' ...']);
    %save(outFile, 'conn_matrix', '-v7.3');


    [conn_matrix, dewhiteM] = icatb_calculate_pca(conn_matrix, sesInfo.userInput.numComp, 'type', 'mpowit', 'whiten', sesInfo.userInput.b_whitening_tmp);
    
    outFile = fullfile(outputDir, conn_dir, [prefix, '_sub', icatb_returnFileIndex(subNum), '_s', num2str(sesNum), '_conn_data.mat']);
    disp(['Saving file ', outFile, ' ...']);
    conn_matrix_ = struct('conn_matrix', conn_matrix);
    save(outFile, '-fromstruct', conn_matrix_);
    disp('Done');
    fprintf('\n');
    
    files(nDataset).name = outFile;
    
end

sesInfo.userInput.numOfSub = numOfSub;
sesInfo.userInput.numOfSess = numOfSess;
sesInfo.userInput.files = files;

%% Save parameter file
paramFile = sesInfo.userInput.param_file;
save(paramFile, 'sesInfo');

