function [beta_weights, eventAvg] = icatb_compute_single_trial_amp(tc, TR, windowSize, selectedOnset, plotResults)
%% Compute single trial amplitudes on fMRI data. This code is an
% implementation of Dr. Eichele's work.
%
% Inputs:
% 1. tc - Scaled ICA timecourse
% 2. TR - Interscan timeinterval
% 3. windowSize - Window size
% 4. selectedOnset - Selected onset
%
% Outputs:
%

icatb_defaults;
global FONT_COLOR;
global UI_FONTUNITS;
global UI_FONTNAME;
global UI_FS;

%% Get event average (Deconvolution method)
[eventAvg, convolved_tc, selectedOnset] = icatb_calculate_eventAvg_deconv(tc, TR, windowSize, selectedOnset);

if (~exist('plotResults', 'var'))
    plotResults = 0;
end

%% Use deconvolution followed by multiple regression to get single trial
% amplitudes. Stimulus onsets are used as regressors (convolved with the
% event average) and convolved timecourse is used as observation.
modelX = zeros(length(tc), length(selectedOnset) + 1);

%% Loop over onsets
for nOns = 1:length(selectedOnset)
    st = zeros(1, length(tc));
    st(selectedOnset(nOns)) = 1;
    st_hrf = conv(st, eventAvg);
    st_hrf = st_hrf(:);
    modelX(:, nOns) = st_hrf(1:length(tc));
end
%% End loop over onsets

modelX(:, end) = 1;

%% Use regression
[beta_weights, R2] = icatb_regress(tc, modelX);

beta_weights = beta_weights(1:end-1);

%% Plot Results
if (plotResults)
    
    figTitle = 'Single trial amplitudes';
    plotTitle = figTitle;
    
    % Figure handle
    graphicsHandle = icatb_getGraphics(figTitle, 'normal', 'single_trial_amplitudes', 'on');
    axesPos = [0.1 0.1 0.8 0.8];
    axesH = axes('parent', graphicsHandle, 'units', 'normalized', 'position', axesPos, 'fontunits', UI_FONTUNITS, 'fontname', ...
        UI_FONTNAME, 'fontSize', UI_FS - 1);
    
    % Plot bar graph
    %bar(tc(selectedOnset), beta_weights(1:end-1), 'parent', axesH);
    plot(selectedOnset, beta_weights, 'parent', axesH);
    
    title(plotTitle, 'color', FONT_COLOR, 'HorizontalAlignment', 'center', 'fontunits', UI_FONTUNITS, 'fontname', UI_FONTNAME, ...
        'fontSize', UI_FS - 1);
    xlabel('Scans', 'parent', axesH);
    ylabel('Amplitudes (\beta)', 'interpreter', 'tex', 'parent', axesH);
    
    set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
    axis(axesH, 'tight');
    
end