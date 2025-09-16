function H = icatb_plot_with_ste_area(H, x, y, e, linecolor, areacolor)
%  Plots a function with error bars
%  If error (e) is not provided,
% USAGE:
%   H = plot_with_ste_area(H, x, y, e, linecolor, areacolor)
% REQUIRED INPUTS:
%   H                   = axis/figure handle
%   x                   = xaxis of plot
%   y                   = data in nreps x points
%   e                   = error to be plotted; if e is [1 x points],
%                           assumes symmetric error, otherwise should be [2 x npoints]
%                           if e is empty, error is the standard error of y, ste = (std(y)/sqrt(nreps))
% OPTIONAL INPUTS:
%   linecolor           = color of line connecting points
%   areacolor           = color of ste area

%
% OUTPUTS:
%   H       = axis handle
%
%  see also: plot_errorbar_area()

if nargin < 3
    y = x;
    x = 1:size(y,2);
end

if nargin < 4
    e = [];
end

if nargin < 5
    linecolor = [0 0 0];
end

if nargin < 6
    areacolor = linecolor + 0.75;
    areacolor(find(areacolor < 0)) = 0;
    areacolor(find(areacolor > 1)) = 1;
end

if isempty(x)
    x = 1:size(y,2);
end

if isvector(y) && length(x)>1
    y = reshape(y, 1, length(y));
end

if length(x) > 1
    goodx = find(~isnan(x));
    x = x(goodx);
    y = y(:,goodx);
end


if isvector(y) && ~isempty(e)
    ymean = y;
elseif isvector(y) && isempty(e) && length(x)>1
    %% x and y are in groups -- need to sort
    y = reshape(y, 1, length(y));
    ux = unique(x);
    ymean = zeros(1,length(ux));
    e = zeros(1, length(ux));
    for ii = 1:length(ux)
        ymean(ii) = mean(y(find(x == ux(ii))));
        e(ii) = std(y(find(x == ux(ii))))/sqrt(length(find(x == ux(ii))));
    end
    x = ux;
end

if isempty(e) | e == 0
    e = std(y)/sqrt(size(y,1));
end

if isempty(H)
    H =   figure;
else
    if strcmp(get(H, 'Type'), 'figure')
        hold on
        figure(H)
    elseif strcmp(get(H, 'Type'), 'axes')
        hold on
        subplot(H)
    end
end

%% begin plotting
if ~exist('ymean', 'var')
    ymean = mean(y,1);
end
icatb_plot_errorbar_area(H, x, ymean, e, linecolor, areacolor);
