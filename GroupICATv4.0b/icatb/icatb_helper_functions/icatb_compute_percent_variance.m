function [ss_residual, ss_total, percent_variance]  = icatb_compute_percent_variance(data, model)
%% Use regression to get explained percent variance in the data by the model
% where components are treated as model.
%
% Inputs:
% 1. data - Data is of size voxels by timepoints
% 2. model - Component timecourses of dimensions timepoints by components
%
% Outputs:
% 1. ss_residual - Sum of squares of residual
% 2. ss_total - Total sum of squares
% 3. percent_variance - Percent variance

%% Initialise variables
voxels = size(data, 1);
time_points = size(data, 2);

%% Do error check on number of timepoints
if (time_points ~= size(model, 1))
    model = model';
end

if (time_points ~= size(model, 1))
    error('Error:DataDimensions', 'Please check the model dimensions as it doesn''t match the data');
end

%% Remove mean of the model
model = icatb_remove_mean(model);

try
    
    %% Arrange data (Timepoints by voxels)
    data = data';
    
    %% Remove mean of the data
    data = icatb_remove_mean(data);
    
    %% Residual
    residual = data - model*(pinv(model)*data);
    
    %% Total sum of squares and residual sum of squares
    ss_total = sum(sum(data.^2), 2);
    ss_residual = sum(sum(residual.^2), 2);
    
catch
    
    %% Initialise vars
    ss_residual = 0;
    ss_total = 0;
    
    %% Loop over in brain voxels
    for nVoxel = 1:voxels    
        
        % Bold signal at current voxel
        bold_signal = data(nVoxel, :)';
        bold_signal = detrend(bold_signal, 0);
        bold_signal = bold_signal(:);    
        
        % Do Regression
        beta_weights = pinv(model)*bold_signal;
        residual = bold_signal - model*beta_weights;    
        
        %% Compute sum of squares and residual sum of squares
        ss_total = ss_total +  sum(bold_signal.^2);
        ss_residual = ss_residual + sum(residual.^2);    
        
    end
    %% End loop over in brain voxels
    
end

%% Percent variance in the data
percent_variance = 100*(1 - (ss_residual/ss_total));
