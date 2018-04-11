function icatb_disp_temp_regress_results(temporal_stats_betas, plotNextPrev)
% Display temporal regression results
%

icatb_defaults;
global FONT_COLOR;
global UI_FS;

if (~exist('plotNextPrev', 'var'))
    plotNextPrev = 1;
end

numOfSess = size(temporal_stats_betas.betas, 1);
numOfSub =  size(temporal_stats_betas.betas, 2);
numComp = length(temporal_stats_betas.component_numbers);
selectedRegressors = temporal_stats_betas.selectedRegressors;
numCond = length(selectedRegressors);
%numOfSess = sesInfo.numOfSess;
ttest_betas = temporal_stats_betas.ttest.tstat;
pval_betas = temporal_stats_betas.ttest.pval;
comps_betas = temporal_stats_betas.component_numbers;
pval_sig1 = 0.01;
pval_sig2 = 0.001;
maxDF = max(squeeze(temporal_stats_betas.ttest.df(:)));

if ((numOfSess == 1) || (numCond == 1))
    figType = 'graphics';
    numFigures = numCond + 1;
    ttest_betas = squeeze(ttest_betas);
    pval_betas = squeeze(pval_betas);
    %numRows = 1;
else
    figType = 'graphics';
    numFigures = (numComp + 1);
    %numRows = 1;
    %numFigures = ceil((numComp + 1)/numRows);
end

fH = repmat(struct('H', []), 1, numFigures);

%for nFig = 1:numFigures

figPos = get(0, 'DefaultFigurePosition');

sz = get(0, 'Screensize');

figPos = [50, 50, sz(3) - 100, max([figPos(4), 0.75*sz(4)])];

%% R^2 summary
% Multiple regression is perfomed using ICA Timecourses as observations and
% SPM design matrix as model timecourses. R^2 values are reported for each
% component and sorted in descending order.
fH(1).H = icatb_getGraphics('Regress values', figType, 'regress_betas', 'on');
set(fH(1).H, 'resize', 'on');
set(fH(1).H, 'position', figPos);
movegui(fH(1).H, 'center');

axesPos = [0.11, 0.2, 0.85, 0.68];
sh = axes('parent', fH(1).H, 'units', 'normalized', 'position', axesPos, 'Xcolor', FONT_COLOR, 'Ycolor', FONT_COLOR, 'Nextplot', 'add', 'fontsize', UI_FS - 2);
bh = bar(temporal_stats_betas.R2(:));
colormap(hsv(1));
set(sh,'XTick',(1:numComp));
set(sh, 'XTickLabel', num2str(temporal_stats_betas.component_numbers(:)));
ylabel('R^2', 'parent', sh, 'interpreter', 'tex', 'color', FONT_COLOR, 'fontsize', UI_FS);
xlabel('Components', 'parent', sh, 'color', FONT_COLOR, 'fontsize', UI_FS);
title('R^2', 'parent', sh, 'interpreter', 'tex', 'color', FONT_COLOR,  'fontsize', UI_FS);
axis(sh, 'tight');
try
    set(sh, 'XTickLabelRotation', 90);
catch
    xticklabel_rotate([],90,[],'interpreter','none', 'Color', FONT_COLOR);
end
%set(sh, 'Xcolor', FONT_COLOR);
%set(sh, 'Ycolor', FONT_COLOR);

if (numOfSub == 1)
    return;
end


tval_sig1 = abs(icatb_spm_invTcdf(pval_sig1/2, maxDF));
tval_sig2 = abs(icatb_spm_invTcdf(pval_sig2/2, maxDF));

YLIMA = [min([ttest_betas(:);-tval_sig2-0.1]), max([ttest_betas(:);tval_sig2 + 0.1])];

%%


if ((numOfSess == 1) || (numCond == 1))
    for nCond = 1:size(ttest_betas, 2)
        
        tmpTT = ttest_betas(:, nCond);
        tmpTP = pval_betas(:, nCond);
        
        nFig = nCond + 1;
        
        %% T-test summary
        % For each run, one sample t-test is computed for each condition
        % and component. T-values are reported. T-values corresponding to p < 0.01 and p < 0.001
        % lines are shown in yellow and green colors.
        fH(nFig).H = icatb_getGraphics('Regress values', figType, 'regress_betas_T', 'on');
        set(fH(nFig).H, 'resize', 'on');
        set(fH(nFig).H, 'position', figPos);
        movegui(fH(nFig).H, 'center');
        
        sh = axes('parent', fH(nFig).H, 'units', 'normalized', 'position', axesPos, 'Xcolor', FONT_COLOR, 'Ycolor', FONT_COLOR, 'Nextplot', 'add',  'fontsize', UI_FS - 2);
        
        bh = bar(tmpTT);
        set(sh, 'Ylim', YLIMA);
        %set(sh,'color', 'w');
        colormap(hsv(1));
        ylabel('T-values', 'parent', sh,  'fontsize', UI_FS);
        xlabel('Components', 'parent', sh,  'fontsize', UI_FS);
        title(selectedRegressors{nCond}, 'parent', sh,  'fontsize', UI_FS);
        set(sh,'XTick',(1:numComp));
        set(sh, 'XTickLabel', num2str(temporal_stats_betas.component_numbers(:)));
        
        hold on;
        xLIMA = get(sh, 'XLim');
        %xLIMA = [1, numComp];
        %set(sh, 'XLim', xLIMA);
        plot(xLIMA, tval_sig1*ones(1, 2), '--y', 'linewidth', 1.5);
        hold on;
        plot(xLIMA, tval_sig2*ones(1, 2), '--g', 'linewidth', 1.5);
        hold on;
        plot(xLIMA, -tval_sig1*ones(1, 2), '--y', 'linewidth', 1.5);
        hold on;
        plot(xLIMA, -tval_sig2*ones(1, 2), '--g', 'linewidth', 1.5);
        hold off;
        %axis(sh, 'tight');
        
        legendStr = {'Comp'; ['p < ', num2str(pval_sig1, '%0.2f')]; ['p < ', num2str(pval_sig2, '%0.3f')]};
        legend(legendStr{:}, 'location', 'best');
        
        try
            set(sh, 'XTickLabelRotation', 90);
        catch
            xticklabel_rotate([],90,[],'interpreter','none', 'Color', FONT_COLOR);
        end
        %end
        
    end
    
else
    
    legendStr = cellstr([repmat('Run ', numOfSess, 1),  num2str((1:numOfSess)')]);
    legendStr = [legendStr; ['p < ', num2str(pval_sig1, '%0.2f')]; ['p < ', num2str(pval_sig2, '%0.3f')]];
    for nComp = 1:numComp
        tmpTP = squeeze(pval_betas(:, nComp, :));
        tmpTT = squeeze(ttest_betas(:, nComp, :));
        nFig = nComp + 1;
        
        fH(nFig).H = icatb_getGraphics('Regress values', figType, 'regress_betas_T', 'on');
        set(fH(nFig).H, 'resize', 'on');
        set(fH(nFig).H, 'position', figPos);
        movegui(fH(nFig).H, 'center');
        
        sh = axes('parent', fH(nFig).H, 'units', 'normalized', 'position', axesPos, 'Xcolor', FONT_COLOR, 'Ycolor', FONT_COLOR, 'Nextplot', 'add',  'fontsize', UI_FS - 2);
        
        bh = bar(tmpTT');
        set(sh, 'Ylim', YLIMA);
        %set(sh,'color', 'w');
        colormap(hsv(numOfSess));
        ylabel('T-values', 'parent', sh,  'fontsize', UI_FS);
        %xlabel('Regressors', 'parent', sh);
        title(['Component ', num2str(temporal_stats_betas.component_numbers(nComp))], 'parent', sh,  'fontsize', UI_FS);
        set(sh,'XTick',(1:length(selectedRegressors)));
        set(sh, 'XTickLabel', selectedRegressors, 'TickDir', 'out');
        
        hold on;
        xLIMA = get(sh, 'XLim');
        %xLIMA = [1, length(selectedRegressors)];
        %set(sh, 'XLim', xLIMA);
        plot(xLIMA, tval_sig1*ones(1, 2), '--y', 'linewidth', 1.5);
        hold on;
        plot(xLIMA, tval_sig2*ones(1, 2), '--g', 'linewidth', 1.5);
        hold on;
        plot(xLIMA, -tval_sig1*ones(1, 2), '--y', 'linewidth', 1.5);
        hold on;
        plot(xLIMA, -tval_sig2*ones(1, 2), '--g', 'linewidth', 1.5);
        hold off;
        try
            set(sh, 'XTickLabelRotation', 90);
        catch
            xticklabel_rotate([],90,[],'interpreter','none', 'Color', FONT_COLOR);
        end
        
        %if (nComp == numComp)
        legend(legendStr{:}, 'location', 'best');
        %end
        %axis(sh, 'tight');
        
    end
    
end
%end


if (plotNextPrev)
    icatb_plotNextPreviousExitButtons(fH);
end


function hText = xticklabel_rotate(XTick,rot,varargin)
%hText = xticklabel_rotate(XTick,rot,XTickLabel,varargin)     Rotate XTickLabel
%
% Syntax: xticklabel_rotate
%
% Input:
% {opt}     XTick       - vector array of XTick positions & values (numeric)
%                           uses current XTick values or XTickLabel cell array by
%                           default (if empty)
% {opt}     rot         - angle of rotation in degrees, 90? by default
% {opt}     XTickLabel  - cell array of label strings
% {opt}     [var]       - "Property-value" pairs passed to text generator
%                           ex: 'interpreter','none'
%                               'Color','m','Fontweight','bold'
%
% Output:   hText       - handle vector to text labels
%
% Example 1:  Rotate existing XTickLabels at their current position by 90?
%    xticklabel_rotate
%
% Example 2:  Rotate existing XTickLabels at their current position by 45? and change
% font size
%    xticklabel_rotate([],45,[],'Fontsize',14)
%
% Example 3:  Set the positions of the XTicks and rotate them 90?
%    figure;  plot([1960:2004],randn(45,1)); xlim([1960 2004]);
%    xticklabel_rotate([1960:2:2004]);
%
% Example 4:  Use text labels at XTick positions rotated 45? without tex interpreter
%    xticklabel_rotate(XTick,45,NameFields,'interpreter','none');
%
% Example 5:  Use text labels rotated 90? at current positions
%    xticklabel_rotate([],90,NameFields);
%
% Example 6:  Multiline labels
%    figure;plot([1:4],[1:4])
%    axis([0.5 4.5 1 4])
%    xticklabel_rotate([1:4],45,{{'aaa' 'AA'};{'bbb' 'AA'};{'ccc' 'BB'};{'ddd' 'BB'}})
%
% Note : you can not RE-RUN xticklabel_rotate on the same graph.
%



% This is a modified version of xticklabel_rotate90 by Denis Gilbert
% Modifications include Text labels (in the form of cell array)
%                       Arbitrary angle rotation
%                       Output of text handles
%                       Resizing of axes and title/xlabel/ylabel positions to maintain same overall size
%                           and keep text on plot
%                           (handles small window resizing after, but not well due to proportional placement with
%                           fixed font size. To fix this would require a serious resize function)
%                       Uses current XTick by default
%                       Uses current XTickLabel is different from XTick values (meaning has been already defined)

% Brian FG Katz
% bfgkatz@hotmail.com
% 23-05-03
% Modified 03-11-06 after user comment
%	Allow for exisiting XTickLabel cell array
% Modified 03-03-2006
%   Allow for labels top located (after user comment)
%   Allow case for single XTickLabelName (after user comment)
%   Reduced the degree of resizing
% Modified 11-jun-2010
%   Response to numerous suggestions on MatlabCentral to improve certain
%   errors.
% Modified 23-sep-2014
%   Allow for mutliline labels


% Other m-files required: cell2mat
% Subfunctions: none
% MAT-files required: none
%
% See also: xticklabel_rotate90, TEXT,  SET

% Based on xticklabel_rotate90
%   Author: Denis Gilbert, Ph.D., physical oceanography
%   Maurice Lamontagne Institute, Dept. of Fisheries and Oceans Canada
%   email: gilbertd@dfo-mpo.gc.ca  Web: http://www.qc.dfo-mpo.gc.ca/iml/
%   February 1998; Last revision: 24-Mar-2003

% check to see if xticklabel_rotate has already been here (no other reason for this to happen)
if isempty(get(gca,'XTickLabel')),
    error('xticklabel_rotate : can not process, either xticklabel_rotate has already been run or XTickLabel field has been erased')  ;
end

% if no XTickLabel AND no XTick are defined use the current XTickLabel
%if nargin < 3 & (~exist('XTick') | isempty(XTick)),
% Modified with forum comment by "Nathan Pust" allow the current text labels to be used and property value pairs to be changed for those labels
if (nargin < 3 || isempty(varargin{1})) & (~exist('XTick') | isempty(XTick)),
    xTickLabels = get(gca,'XTickLabel')  ; % use current XTickLabel
    if ~iscell(xTickLabels)
        % remove trailing spaces if exist (typical with auto generated XTickLabel)
        temp1 = num2cell(xTickLabels,2)         ;
        for loop = 1:length(temp1),
            temp1{loop} = deblank(temp1{loop})  ;
        end
        xTickLabels = temp1                     ;
    end
    varargin = varargin(2:length(varargin));
end

% if no XTick is defined use the current XTick
if (~exist('XTick') | isempty(XTick)),
    XTick = get(gca,'XTick')        ; % use current XTick
end

%Make XTick a column vector
XTick = XTick(:);

if ~exist('xTickLabels'),
    % Define the xtickLabels
    % If XtickLabel is passed as a cell array then use the text
    if (length(varargin)>0) & (iscell(varargin{1})),
        xTickLabels = varargin{1};
        varargin = varargin(2:length(varargin));
    else
        xTickLabels = num2str(XTick);
    end
end

if length(XTick) ~= length(xTickLabels),
    error('xticklabel_rotate : must have same number of elements in "XTick" and "XTickLabel"')  ;
end

%Set the Xtick locations and set XTicklabel to an empty string
set(gca,'XTick',XTick,'XTickLabel','')

if nargin < 2,
    rot = 90 ;
end

% Determine the location of the labels based on the position
% of the xlabel
hxLabel = get(gca,'XLabel');  % Handle to xlabel
xLabelString = get(hxLabel,'String');

% if ~isempty(xLabelString)
%    warning('You may need to manually reset the XLABEL vertical position')
% end

set(hxLabel,'Units','data');
xLabelPosition = get(hxLabel,'Position');
y = xLabelPosition(2);

%CODE below was modified following suggestions from Urs Schwarz
y=repmat(y,size(XTick,1),1);
% retrieve current axis' fontsize
fs = get(gca,'fontsize');

if ~iscell(xTickLabels)
    % Place the new xTickLabels by creating TEXT objects
    hText = text(XTick, y, xTickLabels,'fontsize',fs);
else
    % Place multi-line text approximately where tick labels belong
    for cnt=1:length(XTick),
        hText(cnt) = text(XTick(cnt),y(cnt),xTickLabels{cnt},...
            'VerticalAlignment','top', 'UserData','xtick');
    end
end

% Rotate the text objects by ROT degrees
%set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})
% Modified with modified forum comment by "Korey Y" to deal with labels at top
% Further edits added for axis position
xAxisLocation = get(gca, 'XAxisLocation');
if strcmp(xAxisLocation,'bottom')
    set(hText,'Rotation',rot,'HorizontalAlignment','right',varargin{:})
else
    set(hText,'Rotation',rot,'HorizontalAlignment','left',varargin{:})
end

% Adjust the size of the axis to accomodate for longest label (like if they are text ones)
% This approach keeps the top of the graph at the same place and tries to keep xlabel at the same place
% This approach keeps the right side of the graph at the same place

set(get(gca,'xlabel'),'units','data')           ;
labxorigpos_data = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'ylabel'),'units','data')           ;
labyorigpos_data = get(get(gca,'ylabel'),'position')  ;
set(get(gca,'title'),'units','data')           ;
labtorigpos_data = get(get(gca,'title'),'position')  ;

set(gca,'units','pixel')                        ;
set(hText,'units','pixel')                      ;
set(get(gca,'xlabel'),'units','pixel')          ;
set(get(gca,'ylabel'),'units','pixel')          ;
% set(gca,'units','normalized')                        ;
% set(hText,'units','normalized')                      ;
% set(get(gca,'xlabel'),'units','normalized')          ;
% set(get(gca,'ylabel'),'units','normalized')          ;

origpos = get(gca,'position')                   ;

% textsizes = cell2mat(get(hText,'extent'))       ;
% Modified with forum comment from "Peter Pan" to deal with case when only one XTickLabelName is given.
x = get( hText, 'extent' );
if iscell( x ) == true
    textsizes = cell2mat( x ) ;
else
    textsizes = x;
end

largest =  max(textsizes(:,3))                  ;
longest =  max(textsizes(:,4))                  ;

laborigext = get(get(gca,'xlabel'),'extent')    ;
laborigpos = get(get(gca,'xlabel'),'position')  ;

labyorigext = get(get(gca,'ylabel'),'extent')   ;
labyorigpos = get(get(gca,'ylabel'),'position') ;
leftlabdist = labyorigpos(1) + labyorigext(1)   ;

% assume first entry is the farthest left
leftpos = get(hText(1),'position')              ;
leftext = get(hText(1),'extent')                ;
leftdist = leftpos(1) + leftext(1)              ;
if leftdist > 0,    leftdist = 0 ; end          % only correct for off screen problems

% botdist = origpos(2) + laborigpos(2)            ;
% newpos = [origpos(1)-leftdist longest+botdist origpos(3)+leftdist origpos(4)-longest+origpos(2)-botdist]
%
% Modified to allow for top axis labels and to minimize axis resizing
if strcmp(xAxisLocation,'bottom')
    newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
        origpos(2)+((longest+laborigpos(2))-get(gca,'FontSize')) ...
        origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
        origpos(4)-((longest+laborigpos(2))-get(gca,'FontSize'))]  ;
else
    newpos = [origpos(1)-(min(leftdist,labyorigpos(1)))+labyorigpos(1) ...
        origpos(2) ...
        origpos(3)-(min(leftdist,labyorigpos(1)))+labyorigpos(1)-largest ...
        origpos(4)-(longest)+get(gca,'FontSize')]  ;
end
set(gca,'position',newpos)                      ;

% readjust position of text labels after resize of plot
set(hText,'units','data')                       ;
for loop= 1:length(hText),
    set(hText(loop),'position',[XTick(loop), y(loop)])  ;
end

% adjust position of xlabel and ylabel
laborigpos = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'xlabel'),'position',[laborigpos(1) laborigpos(2)-longest 0])   ;

% switch to data coord and fix it all
set(get(gca,'ylabel'),'units','data')                   ;
set(get(gca,'ylabel'),'position',labyorigpos_data)      ;
set(get(gca,'title'),'position',labtorigpos_data)       ;

set(get(gca,'xlabel'),'units','data')                   ;
labxorigpos_data_new = get(get(gca,'xlabel'),'position')  ;
set(get(gca,'xlabel'),'position',[labxorigpos_data(1) labxorigpos_data_new(2)])   ;


% Reset all units to normalized to allow future resizing
set(get(gca,'xlabel'),'units','normalized')          ;
set(get(gca,'ylabel'),'units','normalized')          ;
set(get(gca,'title'),'units','normalized')          ;
set(hText,'units','normalized')                      ;
set(gca,'units','normalized')                        ;

if nargout < 1,
    clear hText
end
