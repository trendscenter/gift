function [y, points, bandwidth] = icatb_kde(data, points, bandwidth)
%% Kernel density estimator. Using gaussian kernel.
%
% Inputs:
% 1. data - Data must be a vector
% 2. points - Points must be a vector
% 3. bandwith - Bandwidth
%
% Outputs:
% 1. y - density values
% 2. points - Points on which density values are computed
% 3. bandwith - Band width
%

if (numel(data) ~= length(data))
    error('Data must be a vector');
end

% band width
N = length(data);

M = 128;

if (~exist('bandwidth', 'var'))
    med = median(data);
    sigma = median(abs(data-med)) / 0.6745;
    bandwidth = sigma * (4/(3*N))^(1/5);
end

if (~exist('points', 'var'))
    % points
    xmin = min(data) - 3*bandwidth;
    xmax = max(data) + 3*bandwidth;
    points = linspace(xmin, xmax, M);
end

if (numel(points) ~= length(points))
    error('Data must be a vector');
end

M = length(points);

%% Kernel density estimator
if (~isempty(which('ksdensity.m')))
    [y, points] = ksdensity(data, points, 'width', bandwidth);
else
    if (~isempty(which('bsxfun.m')))
        y = mean(icatb_norm_pdf((1/bandwidth) * bsxfun(@minus, points, data(:)))) ./ bandwidth;
    else
        y = mean(icatb_norm_pdf((1/bandwidth) * (repmat(points(:)', N, 1) - repmat(data(:), 1, M)))) ./ bandwidth;
    end
end