function [refData] = icatb_construct_refData(selectedRegressors, numDataSets, numRegress)
% construct refData structure with reference function names

for ii = 1:numDataSets
    startTp = numRegress*(ii - 1) + 1;
    endTp = numRegress*ii;
    refNames = deblank(selectedRegressors(startTp:endTp, :));
    refData(ii).refNames = refNames;
    clear refNames;
end