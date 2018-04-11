function [estimateComp, estCompVec, meanMdl] = icatb_parComputeMDL_cluster(sesInfo, mask, precisionType)
%% Compute MDL in parallel mode

estimateComp = 0;

disp('Computing MDL for all subjects and sessions in parallel. Please wait ...');

% Load functional data files here and mask out voxels
% and then call estimate dimension to return the number of components
estCompVec = zeros( length(sesInfo.userInput.files), 1 );
tpMin = min(sesInfo.userInput.diffTimePoints) - 1;
meanMdl = zeros(1, tpMin);
files = sesInfo.userInput.files;
drawnow;

parfor numFiles = 1:length(sesInfo.userInput.files)

    fprintf('\n');
    disp('..............................................');
    disp('');
    helpHandleName = ['Loading data set ', num2str(numFiles)];
    disp(helpHandleName);
    drawnow;
    [ estCompVec(numFiles), mdlVal] = icatb_estimate_dimension(files(numFiles).name, mask, precisionType);

    dd = min([tpMin, length(mdlVal)]);
    tmpM = zeros(1, tpMin);
    tmpM(1:dd) = mdlVal(1:dd);

    meanMdl = meanMdl + tmpM;

    estimateComp = estCompVec(numFiles) + estimateComp;

    drawnow;

end

meanMdl = meanMdl/length(sesInfo.userInput.files);

estimateComp = round(estimateComp/length(sesInfo.userInput.files));

