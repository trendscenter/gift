function [tc, spatial_maps] = icatb_dual_regress(y, X)
%% Use dual regression approach to compute the spatial maps and time
% courses
%
% Inputs:
% 1. y - Observations in columns (Voxels by time points)
% 2. X - design matrix (Voxels by components)
%
% Outputs:
% 1. tc - Time courses (Timepoints by components)
% 2. spatial_maps - Spatial maps (Components by voxels)
%

%% First step. Fit model matrix to the data to get time courses.
X = icatb_remove_mean(X);

tc = pinv(X)*icatb_remove_mean(y);
tc = tc';
clear X;

% Store mean of timecourses
mean_tc = mean(tc);

% Remove mean of timecourse
tc = icatb_remove_mean(tc);

%% Second step. Fit Time courses at each voxel to get the spatial maps.
try

    spatial_maps = pinv(tc)*icatb_remove_mean(y');

catch
    %% Use less memory usage way to do regression

    % Initialise spatial maps
    spatial_maps = zeros(size(tc, 2), size(y, 1));


    % Loop over voxels
    for nVoxel = 1:size(y, 1)

        bold_signal = detrend(y(nVoxel, :), 0);

        spatial_maps(:, nVoxel) = pinv(tc)*bold_signal(:);

        clear bold_signal;

    end
    % End of loop over voxels

end

clear y;

%% Add mean back to the timecourses
tc = tc + (repmat(mean_tc, size(tc, 1), 1));
