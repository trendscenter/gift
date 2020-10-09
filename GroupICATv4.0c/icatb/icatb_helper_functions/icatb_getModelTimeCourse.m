function [modelX, model_inds] = icatb_getModelTimeCourse(tp, spmInfo, selectedRegressors, nSess)
%% GET SPM model timecourse

modelX = NaN(tp, length(selectedRegressors));


XX = spmInfo.SPM.xX.X;
names = cellstr(char(spmInfo.SPM.xX.name));

try
    spmInfo.SPM.nscan(nSess);
catch
    nSess = 1;
end

if (nSess == 1)
    startTp = 1;
else
    startTp = sum(spmInfo.SPM.nscan(1:nSess-1)) + 1;
end

endTp = sum(spmInfo.SPM.nscan(1:nSess));
model_inds = [];
for nR = 1:length(selectedRegressors)
    regressName = ['Sn(', num2str(nSess), ') ', selectedRegressors{nR}];
    inds = strmatch(lower(regressName), lower(spmInfo.SPM.xX.name), 'exact');
    if (~isempty(inds))
        model_inds = [model_inds, nR];
        modelX(:, nR) = XX(startTp:endTp, inds(1));
    end
end