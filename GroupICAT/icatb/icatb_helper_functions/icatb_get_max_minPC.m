function [maxPC, minPC] = icatb_get_max_min(diffTimePoints, numRedSteps, PC1, PC2, PC3)

if min(diffTimePoints) < 2
    error('Need atleast 2 images to proceed through group ICA');
end

if min(PC1) < 2 | min(PC2) < 2 | min(PC3) < 2
    error('Need atleast two principal components');
end

if numRedSteps == 1
    maxPC = min(diffTimePoints);
    minPC = 2;

elseif numRedSteps == 2
    maxPC = 4*PC1;
    minPC = 2;
else
    maxPC = ceil(length(diffTimePoints) / 4)*PC2;
    minPC = 2;
end