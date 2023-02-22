function H = icatb_simtb_pcolor(X,Y,IM)
%   simtb_pcolor()  - Produces pcolor-like images of data matrices
%
%   Usage:
%    >> H = simtb_pcolor(X, Y, IM);
%
%   INPUTS: 
%   X            = [1 x m] vector of labels (ticks) in the x-dimension
%   Y            = [1 x n] vector of labels (ticks) in the y-dimension 
%   IM           = [m x n] matrix of data
%
%   OUTPUTS:
%   H            = Image handle 
%
%   see also: simtb_figure_params

%% using imagesc()
H = imagesc(X,Y,IM);
set(gca, 'ticklength', [0 0])
axis xy
hold on;
xlim = get(gca, 'xlim');
ylim = get(gca, 'ylim');
xx = [xlim(1):xlim(2)]; % where to draw grid lines x 
yy = [ylim(1):ylim(2)]; %   and y
plot(repmat(xx,2,1), repmat(ylim',1,length(xx)), 'k'); % vertical lines
plot(repmat(xlim',1,length(yy)), repmat(yy,2,1), 'k'); % horizontal lines

set(gca,'XTick',X(1:end));
set(gca,'XTickLabel',X(1:end));
set(gca,'YTick',Y(1:end));
set(gca,'YTickLabel',Y(1:end));
