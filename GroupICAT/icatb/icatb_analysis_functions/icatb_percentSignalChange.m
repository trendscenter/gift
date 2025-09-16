function [icasig, A] = icatb_percentSignalChange(origData, icasig, A)
% scales to percent signal change

icatb_defaults;
global DETRENDNUMBER;

% detrend number
detrendNumber = DETRENDNUMBER;

A = A';

%scale components and timecourses
for ii=1:size(icasig,1)
    Bweighted = 0;
    voxelValueTotal = 0;
    % scale using weighted average
    numberOfTopVoxels = 5;
    [sortedICASig sortedIndex] = sort(abs(icasig(ii,:)));
    topFenceStart = size(icasig, 2) - numberOfTopVoxels;
    topVoxels  = sortedICASig([topFenceStart+1:end]);
    topIndex = sortedIndex([topFenceStart+1:end]);
    % detrend time course once
    A(:, ii) = (icatb_detrend(A(:, ii), 1, length(A(:, ii)), detrendNumber)); 
    for j = 1:numberOfTopVoxels;
        %voxel's original BOLD signal
        origVox = origData(topIndex(j), :);           
        origVoxScaled = 100*origVox / mean(origVox);
        origVoxScaled = (icatb_detrend(origVoxScaled, 1, length(origVoxScaled), detrendNumber));      
        % calculate beta coefficients
        [B] = icatb_regress(origVoxScaled, A(:, ii));               
        if isnan(B(1))
            error('Cannot find voxels in the brain. Please specify a suitable mask');
        end
        %weighted regression values
        Bweighted = ((B(1)) * topVoxels(j)) + Bweighted;
        voxelValueTotal = topVoxels(j) + voxelValueTotal;
        clear temp;
    end
    Bweighted = Bweighted/voxelValueTotal;
    %rescale images...(don't change the sign)
    icasig(ii,:) = icasig(ii,:)*(abs(Bweighted/topVoxels(j)));
    %weight tc to maximum percent signal change...(don't change the sign)
    A(:, ii) = (abs(Bweighted)*A(:, ii));
end
