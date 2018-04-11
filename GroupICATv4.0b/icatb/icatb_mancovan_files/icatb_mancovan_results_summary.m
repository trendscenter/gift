function icatb_mancovan_results_summary(mancovanInfo)


try
    structFile = mancovanInfo.structFile;
catch
    structFile = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');
end

try
    imStr = mancovanInfo.image_values;
catch
    imStr = 'Positive';
end

try
    threshdesc = mancovanInfo.threshdesc;
catch
    threshdesc = 'fdr';
end

mancovanInfo.threshdesc = threshdesc;

mancovanInfo.image_values = imStr;
mancovanInfo.structFile = structFile;

%% Features display (Spatial maps, Timecourses spectra, FNC correlations)
% * *a) Spatial maps and spectra* - Orthogonal slices of T-maps for each component and timecourses spectra is shown.
% * *b) FNC correlations* -  FNC correlations are averaged across subjects.
%
gH = icatb_plot_mancova_features(mancovanInfo);

drawnow;

%% Multivariate results
%  Multivariate tests are done on the features to determine the significant covariates which are later used in the univariate tests on each feature.
designCriteria = 'mancova';
try
    designCriteria = mancovanInfo.designCriteria;
catch
end

if (isfield(mancovanInfo, 'univInfo') && ~isempty(mancovanInfo.univInfo))
    designCriteria = 'none';
end

if (strcmpi(designCriteria, 'mancova'))
    gH(1).H = icatb_plot_mult_mancovan(mancovanInfo);
    if (isfield(mancovanInfo, 'time'))
        % Time 1
        gH(end + 1).H = icatb_plot_mult_mancovan(mancovanInfo, 1);
        % Time 2
        gH(end + 1).H  = icatb_plot_mult_mancovan(mancovanInfo, 2);
    end
end

drawnow;

%% Univariate results
%
% * *a) Spatial maps* - T-maps of the significant covariate are shown as composite t-maps. Beta-values are averaged over significant clusters.
% * *b) Spectra* - Univariate t-tests are done using the significant covariates on the spectra. Beta-values ate averaged over frequency bands.
% * *c) FNC* - Univariate t-tests are done using the significant covariates on the FNC correlations.
% Connectogram of FNC correlations is also shown. Thumbnails of mean component maps are also plotted.
outDir = mancovanInfo.outputDir;
outputFiles = mancovanInfo.outputFiles;

start_terms = {};
for nO = 1:length(outputFiles)
    for nR = 1:length(outputFiles(nO).filesInfo.result_files)
        load(fullfile(outDir, outputFiles(nO).filesInfo.result_files{nR}), 'UNI');
        for nT = 1:length(UNI.tests)
            for nC = 1:length(UNI.stats{nT}.Contrast)
                if (~isempty(UNI.stats{nT}.Contrast{nC}))
                    start_terms{end + 1} = UNI.stats{nT}.Contrast{nC};
                else
                    start_terms{end + 1} = UNI.tests{nT};
                end
            end
        end
    end
end

if (~isempty(start_terms))
    [ddd, indbb] = unique(start_terms);
    start_terms = start_terms(sort(indbb));
end

load(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(1).filesInfo.result_files{1}), 'UNI');
if (~exist('UNI', 'var') || isempty(UNI))
    error('Please run mancova in order to view results');
end

for nF = 1:length(start_terms)
    mancovanInfo.covariatesToPlot = start_terms(nF);
    gH  = icatb_plot_univariate_results(mancovanInfo);
    drawnow;
    
    if (isfield(mancovanInfo, 'time'))
        % Time 1
        gH = icatb_plot_univariate_results(mancovanInfo, 1);
        drawnow;
        % Time 2
        gH  =  icatb_plot_univariate_results(mancovanInfo, 2);
        drawnow;
    end
end


function applyDefs(figH)

for nH = 1:length(figH)
    
    H = figH(nH).H;
    
    pos = get(0, 'defaultFigurePosition');
    set(H, 'color', 'w');
    %set(H, 'position', pos);
    set(H, 'resize', 'on');
    set(H, 'visible', 'on');
    titleH = get(findobj(H, 'type', 'axes'), 'title');
    if (iscell(titleH))
        % titleH = cell2mat(titleH);
        for nT = 1:length(titleH)
            set(titleH{nT}, 'color', 'k');
        end
    else
        set(titleH, 'color', 'k');
    end
    axesH = findobj(H, 'type', 'axes');
    set(axesH, 'YColor', 'k', 'XColor', 'k');
    set(axesH, 'fontname', 'times');
    set(axesH, 'fontsize', 12);
    textH = findobj(H, 'type', 'text');
    set(textH, 'color', 'k');
    set(textH, 'fontname', 'times');
    set(textH, 'fontsize', 12);
    try
        C = findobj(H, 'type', 'colorbar');
        set(C, 'color', 'k');
        set(get(C, 'label'), 'color', 'k');
    catch
    end
    
end
