function gH = icatb_plot_mult_mancovan(mancovanInfo, timeNo)
%% Plot Multi variate results
%
% Inputs:
% 1. mancovanInfo - Mancovan information file
%

if (~exist('timeNo', 'var'))
    timeNo = 0;
end

outputDir = mancovanInfo.outputDir;
comps = mancovanInfo.comps;
res = mancovanInfo.outputFiles;

icatb_defaults;
global UI_FONTNAME;
global UI_FS;

load icatb_colors redV;

% fig_title = 'Multivariate Results';
% gH = icatb_getGraphics(fig_title, 'graphics', 'mult_variate_results', 'on');
% nrows = ceil(length(res)/2);
% ncols = ceil(length(res)/nrows);

%load icatb_colors hot;

%gH = zeros(1, length(res));

for nRes = 1:length(res)
    %sh = subplot(nrows, ncols, nRes);
    %results = res{nRes};
    
    fig_title = ['Multivariate Results (', res(nRes).feature_name, ')'];
    gH(nRes) = icatb_getGraphics(fig_title, 'graphics', 'mult_variate_results', 'on');
    colormap(jet(64));
    sh = axes('parent', gH(nRes), 'units', 'normalized', 'position', [0.16, 0.16, 0.74, 0.74]);
    
    result_files = res(nRes).filesInfo.result_files;
    %comps = [];
    for nR = 1:length(result_files)
        if (timeNo > 0)
            load(fullfile(outputDir, result_files{nR}), 'time');
            MULT = time.MULT{timeNo};
        else
            load(fullfile(outputDir, result_files{nR}), 'MULT');
        end
        final_terms = MULT.final_terms;
        pvalues = MULT.p;
        start_terms = MULT.start_terms;
        
        if (nR == 1)
            STATS = NaN*ones(length(result_files), length(start_terms));
        end
        
        for nT = 1:length(final_terms)
            index = strmatch(final_terms{nT}, start_terms, 'exact');
            STATS(nR, index) = pvalues(nT);
        end
        
        clear MULT;
        
    end
    
    
    %% Plot output
    featureName = upper(res(nRes).feature_name);
    
    if (strcmpi(featureName, 'fnc correlations'))
        comp_str = {'ALL'};
        if (length(result_files) == 2)
            comp_str = {'ALL Comps'; 'Networks Avg'};
        end
    elseif (strcmpi(featureName, 'fnc correlations (lag)'))
        comp_str = {'Correlations'; 'Lag'};
    else
        try
            comp_str = cellstr(num2str(comps(:)))';
        catch
        end
    end
    
    logp = -log10(STATS + eps);
    
    cmap = redV; %cmap = cmap(1:size(hot, 1)*0.5,:); cmap = sqrt(cmap);
    
    xlabelStr = 'Components';
    ylabelStr = '-log_1_0(p-value)';
    ytickoff = 0;
    %     if (mod(nRes, 2) ~= 0)
    %         ytickoff = 0;
    %     end
    
    if (timeNo > 0)
        titleStr = ['Mancovan(Time', num2str(timeNo), ')_', featureName];
        titleName = [featureName, '(Time ', num2str(timeNo), ')'];
    else
        titleStr = ['Mancovan_', featureName];
        titleName = featureName;
    end
    
    vars = {logp, comp_str, start_terms, 'fig_title', titleStr, 'tag', titleStr, 'title', titleName, 'cmap', cmap, 'clim',  [0 -log10(eps)], 'xlabel', ...
        xlabelStr, 'ylabel', ylabelStr, 'axesH', sh, 'ytickoff', ytickoff};
    [dd, ch] = icatb_plot_matrix(vars{:});
    
    childAxesH = get([sh, ch], 'children');
    for nChildH = 1:length(childAxesH)
        if iscell(childAxesH)
            set(childAxesH{nChildH}, 'hittest', 'off');
        else
            set(childAxesH(nChildH), 'hittest', 'off');
        end
    end
    %     if (iscell(childAxesH))
    %         childAxesH = cell2mat(childAxesH);
    %     end
    %     set(childAxesH, 'hittest', 'off');
    vars{end} = 0;
    %set([sh, ch], 'ButtonDownFcn', {@openMatrix, vars});
    clear vars;
    
end

%
% axisH = axes('Parent', gH, 'position', [0 0 1 1], 'visible', 'off');
% xPos = 0.5; yPos = 0.97;
% titleColor = 'c';
% text(xPos, yPos, fig_title, 'color', titleColor, 'fontweight', 'bold', 'fontsize', UI_FS + 2, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME, 'parent', axisH);


function openMatrix(hObject, event_data, vars)
%% Open Matrix
%

figH = get(hObject, 'parent');

if (strcmpi(get(figH, 'selectionType'), 'open'))
    set(figH, 'pointer', 'watch');
    gH = icatb_getGraphics(get(figH, 'name'), 'normal', 'Mult results', 'off');
    colormap(get(figH, 'colormap'));
    sh = subplot(1, 1, 1);
    for nV = 1:length(vars)
        if (ishandle(vars{nV}))
            break;
        end
    end
    vars{nV} = sh;
    icatb_plot_matrix(vars{:});
    set(gH, 'visible', 'on');
    set(gH, 'resize', 'on');
    set(figH, 'pointer', 'arrow');
end