function H = icatb_plot_errorbar_area(H, x, ymean, e, linecolor, areacolor)

if nargin < 5
    linecolor = [0 0 0];
    areacolor = [.5 .5 .5];
end


%% main data
if ~isempty(linecolor)
    K = plot(x,ymean);
    set(K, 'Color', linecolor, 'LineStyle', '-', 'LineWidth', 2);
end

hold on

% axis tight
% if strcmp(get(H, 'Type'), 'figure')
%     drawaxismargins(get(H, 'CurrentAxes'));
% else
%     drawaxismargins(H);
% end



%% error bar may not be symmetric

if size(e,1) > 1
    eup = e(1,:);
    edown = e(2,:);
else
    eup = e(1,:);
    edown = e(1,:);
end

xplot = [x, fliplr(x)];
yplot = [eup+ymean, fliplr(ymean-edown)];
F = fill(xplot,yplot, areacolor);

%% uncomment here to make error bar areas blend
set(F, 'LineStyle', 'none', 'FaceAlpha', 0.5);


% if strcmp(get(H, 'Type'), 'figure')
%     drawaxismargins(get(H, 'CurrentAxes'));
% else
%     drawaxismargins(H);
% end

if ~isempty(linecolor)
    uistack(K, 'top')
end
