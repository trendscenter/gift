function icatb_parComputeMDL(param_file, dataSetsToRun, precisionType, outName)
%% Compute MDL in parallel mode
%

%% Initialise parameters
load (param_file);
estCompVec = zeros(1, length(dataSetsToRun));
tpMin = min(sesInfo.userInput.diffTimePoints) - 1;
meanMdl = zeros(1, tpMin);
mask = zeros(sesInfo.userInput.HInfo.DIM(1:3));
mask(sesInfo.userInput.mask_ind) = 1;
mask = (mask == 1);

files = sesInfo.userInput.files;

%% Run MDL on selected data-sets
count = 0;
for numFiles = dataSetsToRun
    
    count = count + 1;
    fprintf('\n');
    disp('..............................................');
    disp('');
    helpHandleName = ['Loading data set ', num2str(numFiles)];
    disp(helpHandleName);
    drawnow;
    [estCompVec(count), mdlVal] = icatb_estimate_dimension(files(numFiles).name, mask, precisionType);
    meanMdl = meanMdl + mdlVal(1:tpMin);
    drawnow;
    
end

%% Save info
icatb_save(outName, 'estCompVec', 'meanMdl');

