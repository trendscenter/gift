function icatb_setup_dfc_roi
%% Setup dynamic FC for ROI-ROI or ROI-voxel
%

% if (isempty(which('spm')))
%     error('SPM is required for ROI dFC to run. Please add SPM (SPM5 and above) on path');
% end

dfcRoiInfo = setup_dyn_fc;

%% Parse files list and get the timepoints count
filesList = cellstr(dfcRoiInfo.userInput.filesInfo.filesList);
time_points = zeros(1, length(filesList));
for n = 1:length(time_points)
    
    
    currentFileN = deblank(filesList{n});
    [pathstr, fN, extn] = fileparts(currentFileN);
    filesP = icatb_listFiles_inDir(pathstr, [fN, extn]);
    if (isempty(filesP))
        error(['Files doesn''t exist. Please check the file pattern ', currentFileN]);
    end
    filesP = icatb_fullFile('directory', pathstr, 'files', filesP);
    filesP = icatb_rename_4d_file(filesP);
    
    tp = (1:size(filesP, 1));
    if ~isempty(dfcRoiInfo.userInput.filesInfo.file_numbers)
        file_nums = dfcRoiInfo.userInput.filesInfo.file_numbers;
        file_nums(file_nums > max(tp)) = [];
        tp = tp(file_nums);
    end
    
    time_points(n) = length(tp);
    
    filesList{n} = deblank(filesP(tp, :));
    
end

dfcRoiInfo.userInput.filesInfo.filesList = filesList;
dfcRoiInfo.userInput.time_points = time_points;


%% Coregister atlas and masks files to the functional images
refFile = deblank(dfcRoiInfo.userInput.filesInfo.filesList{1}(1,:));
refFile = icatb_rename_4d_file(refFile);
refFile = deblank(refFile(1,:));
noisecloud_spm_coregister(refFile, dfcRoiInfo.userInput.atlas_mask_file, '', dfcRoiInfo.userInput.outputDir);

[pa, fN, extn] = fileparts(dfcRoiInfo.userInput.atlas_mask_file);
dfcRoiInfo.userInput.atlas_mask_file = fullfile(dfcRoiInfo.userInput.outputDir, ['r', fN, extn]);
if (~isempty(dfcRoiInfo.userInput.maskFile))
    VR = icatb_spm_vol(refFile);
    VM = icatb_spm_vol(dfcRoiInfo.userInput.maskFile);
    if (all(VR(1).dim(1:3) ~= VM(1).dim(1:3)))
        [pa, fN, extn] = fileparts(dfcRoiInfo.userInput.maskFile);
        noisecloud_spm_coregister(refFile, dfcRoiInfo.userInput.maskFile, '', dfcRoiInfo.userInput.outputDir);
        dfcRoiInfo.userInput.maskFile = fullfile(dfcRoiInfo.userInput.outputDir, ['r', fN, extn]);
    end
end

V = icatb_spm_vol(refFile);
V.n(1)=1;
dfcRoiInfo.userInput.V = V;

% Save results
disp('Saving parameters ...');
param_file = [dfcRoiInfo.userInput.prefix, '_dfc_roi.mat'];
dfcRoiInfo.userInput.param_file = param_file;
param_file = fullfile(dfcRoiInfo.userInput.outputDir, param_file);
save(param_file, 'dfcRoiInfo');
disp(['Parameters saved in file ', param_file]);
disp('');

%% Run dFC ROI
icatb_run_dfc_roi(dfcRoiInfo);
