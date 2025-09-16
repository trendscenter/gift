function icatb_run_spatial_dfnc(sdfncInfo)
%% Run spatial dfnc
%

%% Initialise variables
outputDir = sdfncInfo.userInput.outputDir;
prefix = sdfncInfo.userInput.prefix;
numComp = sdfncInfo.userInput.numComp;
window_length = sdfncInfo.userInput.window_size;
timepoints = sdfncInfo.userInput.timepoints;

files = sdfncInfo.userInput.files;
subjects = [sdfncInfo.userInput.group.val];
numOfSub = length(subjects);
numOfSess = sdfncInfo.userInput.numOfSess;
numGroups = length(sdfncInfo.userInput.group);
windows = ceil(timepoints/window_length);
windows = windows + (windows - 1);
increments = ceil(window_length/2);

mask_ind = sdfncInfo.userInput.mask_ind;
voxels = length(mask_ind);
sdfncInfo.subjects = subjects;
sdfncInfo.mask_ind = mask_ind;
sdfncInfo.numOfSub = numOfSub;
sdfncInfo.numGroups = numGroups;
sdfncInfo.windows = windows;
sdfncInfo.numComp = numComp;
sdfncInfo.outputDir = outputDir;
sdfncInfo.numOfSess = numOfSess;

num_iva_runs = 1;
try
    num_iva_runs = sdfncInfo.userInput.num_iva_runs;
catch
end

IVA_Options = {};
try
    IVA_Options = sdfncInfo.userInput.IVA_Options;
catch
end

drawnow;


tic;
logFile = fullfile(outputDir, [prefix, '_results.log']);
diary(logFile);

disp('----------------------------------------------------');
disp('-------------- STARTING SPATIAL DYNAMIC FNC --------------');
disp('----------------------------------------------------');


icc = zeros(numComp, voxels, numOfSub*numOfSess*windows);

countNW = 0;
%% Loop over subjects
for nSub = 1:numOfSub
    
    %% Loop over sessions
    for nSess = 1:numOfSess
        
        datasetNo = (subjects(nSub) - 1)*numOfSess + nSess;
        disp(['Loading subject ', num2str(subjects(nSub)), ' session ', num2str(nSess), ' ...']);
        data = icatb_remove_mean(icatb_read_data(files(datasetNo).name, [], mask_ind));
        
        dewhiteM = cell(1, windows);
        
        %% Loop over windows
        for nwin = 1:windows
            
            disp(['Computing PCA on subject ', num2str(subjects(nSub)), ' session ', num2str(subjects(nSess)), ' window ', num2str(nwin), ' ...']);
            
            countNW = countNW + 1;
            
            s = increments*(nwin - 1) + 1;
            e = s + window_length - 1;
            if (nwin == windows)
                e = timepoints;
            end
            
            [tmp, dewhiteM{nwin}] = icatb_calculate_pca(data(:, s:e), numComp);
            icc(:, :, countNW) = tmp';
            fprintf('\n');
        end
        %% End of loop over windows
        
        clear data;
        
        fname = [prefix, '_sdfnc_sub_', icatb_returnFileIndex(subjects(nSub)), '_sess_', icatb_returnFileIndex(nSess), '_results.mat'];
        icatb_save(fullfile(outputDir, fname), 'dewhiteM');
        
        fprintf('\n');
        drawnow;
        
    end
    %% End of loop over sessions
    
    fprintf('\n');
    
end
%% End of loop over subjects

%% DO IVA and backreconstruct
disp('Computing IVA-GL on windowed data-sets ...');
[W, icc] = runIVA(icc, num_iva_runs, IVA_Options);

% [~, W] = icatb_icaAlgorithm('iva-gl', icc);
% for i = 1:size(icc, 3)
%     icc(:, :, i) = squeeze(W(:, :, i))*icc(:, :, i);
% end
disp('Done');
fprintf('\n');

icc = reshape(icc, numComp, voxels, numOfSub, numOfSess, windows);
W = reshape(W, numComp, numComp, numOfSub, numOfSess, windows);

%% Write output
disp('Writing IVA-GL components ...');

outputFiles = cell(numOfSub, numOfSess);
% Loop over subjects
for nSub = 1:numOfSub
    % Loop over sessions
    for nSess = 1:numOfSess
        fname = [prefix, '_sdfnc_sub_', icatb_returnFileIndex(subjects(nSub)), '_sess_', icatb_returnFileIndex(nSess), '_results.mat'];
        ic = squeeze(icc(:, :, nSub, nSess, :));
        weights = squeeze(W(:, :, nSub, nSess, :));
        outputFiles{nSub, nSess} = fname;
        icatb_save(fullfile(outputDir, fname), 'ic', 'weights', '-append');
        clear ic weights;
    end
    % End of loop over sessions
end
% End of loop over subjects

%% End of writing output images

sdfncInfo.resultFiles = outputFiles;
clear outputFiles;


%% Write t-maps across all subjects, sessions and windows
icc = reshape(icc, sdfncInfo.numComp, length(sdfncInfo.mask_ind),  numOfSub*numOfSess*windows);
icc = permute(icc, [3, 2, 1]);

[V, HInfo] = icatb_returnHInfo(deblank(sdfncInfo.userInput.files(1).name(1,:)));
V = V(1);
V.dt(1) = 16;
sdfncInfo.HInfo = HInfo;
fname = [sdfncInfo.userInput.prefix, '_sdfnc_tmap_all.nii'];
V.fname = fullfile(sdfncInfo.userInput.outputDir, fname);
V.n(1) = 1;


for nComp = 1:size(icc, 3)
    V.n(1) = nComp;
    dat = squeeze(icc(:, :, nComp));
    [t, p, stats_u] = mT(dat, ones(size(dat, 1), 1), [], 0, {'verbose'});
    tmp = zeros(V(1).dim(1:3));
    tmp(sdfncInfo.mask_ind) = t;
    icatb_write_vol(V, tmp);
end

sdfncInfo.compFiles = fname;

outputFile = fullfile(outputDir, [prefix, '_sdfnc.mat']);
icatb_save(outputFile, 'sdfncInfo');
fprintf('\n');
disp(['Spatial dFNC information is saved in file ', outputFile]);

fprintf('\n');

totalTime = toc;

disp(['Total time taken to complete the analysis is ', num2str(totalTime/60), ' minutes']);

fprintf('\n');

disp('----------------------------------------------------');
disp('-------------- Finished spatial dFNC Analysis --------------');
disp('----------------------------------------------------');

diary('off');


function [W, icasig] = runIVA(data, num_iva_runs, IVA_Options)

if (num_iva_runs > 1)
    % MST
    icasigR = cell(1, num_iva_runs);
    fprintf('\n');
    disp(['Number of times ICA will run is ', num2str(num_iva_runs)]);
    parfor nRun = 1:length(icasigR)
        fprintf('\n');
        disp(['Run ', num2str(nRun), ' / ', num2str(num_iva_runs)]);
        fprintf('\n');
        [dd1, W, dd3, icasigR{nRun}]  = icatb_icaAlgorithm('iva-gl', data, IVA_Options);
    end
    clear dd1 dd2 dd3;
    
    [corrMetric, W, A, icasig, bestRun] = icatb_bestRunSelection(icasigR, data);
    
else
    [dd1, W, A, icasig]  = icatb_icaAlgorithm('iva-gl', data, IVA_Options);
end

