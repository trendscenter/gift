function Correlations = icatb_partial_corr(currentTC, modelX)
%% Compute partial correlations
%
% Inputs:
% 1. currentTC - Timecourses
% 2. modelX - Regressors of interest
%
% Outputs:
% Correlations - Correlations of dimensions number of regressors by number
% of columns of timecourses

currentTC = icatb_remove_mean(currentTC);
modelX = icatb_remove_mean(modelX);

betas = icatb_regress(currentTC, modelX);
Correlations = zeros(size(modelX, 2), size(currentTC, 2));
for nC = 1:size(modelX, 2)
    comp_inds = (1:size(modelX, 2));
    other_inds = find(comp_inds ~= nC);
    tempTC = currentTC - modelX(:, other_inds)*betas(other_inds, :);
    [partialCorrSlopes, partial_rsquare] = icatb_regress(tempTC, modelX(:, nC));
    Correlations(nC, :) = sign(partialCorrSlopes).*sqrt(abs(partial_rsquare));
end