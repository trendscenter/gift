function [data, XYZ] = icatb_read_vols(V)
% Read data using SPM8 functions

% Initialise data
data = zeros([V(1).dim(1:3), length(V)]);

% read the data
for ii = 1:length(V)
    [data(:, :, :, ii), XYZ] = icatb_spm_read_vols(V(ii));
end

