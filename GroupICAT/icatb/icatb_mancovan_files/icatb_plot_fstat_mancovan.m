function figH = icatb_plot_fstat_mancovan(mancovanInfo)
%% Generate display for F-stat
%

outputDir = mancovanInfo.outputDir;
fname = fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(1).filesInfo.result_files{1});
load(fname);
if (~exist('UNI', 'var') || isempty(UNI))
    return
end

if (~isfield(UNI, 'FstatInfo'))
    figH = [];
    return
end

try
    p_threshold = mancovanInfo.display.p_threshold;
catch
    p_threshold = mancovanInfo.userInput.p_threshold;
end

try
    mancovanInfo.display.t_threshold;
catch
    mancovanInfo.display.t_threshold = 1.5;
end

mancovanInfo.display.p_threshold = p_threshold;

display_connectogram = 1;
try
    display_connectogram = mancovanInfo.display.display_connectogram;
catch
end

mancovanInfo.display.display_connectogram = display_connectogram;

prefix = mancovanInfo.prefix;
threshdesc = 'none';
try
    threshdesc = mancovanInfo.display.threshdesc;
catch
end
mancovanInfo.display.threshdesc = threshdesc;

structFile = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet_3x3x3.nii');
try
    structFile = mancovanInfo.display.structFile;
catch
end

regressors = UNI.FstatInfo.regressors;

image_values = 'positive';
try
    image_values = mancovanInfo.display.image_values;
catch
end

load icatb_colors coldhot;


comp_network_names = cell(length(mancovanInfo.comp), 2);
for nC = 1:size(comp_network_names, 1)
    comp_network_names{nC, 1} = mancovanInfo.comp(nC).name;
    comp_network_names{nC, 2} = mancovanInfo.comp(nC).value;
end

figH = [];
for nFeature = 1:length(mancovanInfo.outputFiles)
    
    feature_name = capitalize(mancovanInfo.outputFiles(nFeature).feature_name);
    result_files = mancovanInfo.outputFiles(nFeature).filesInfo.result_files;
    
    if (strcmpi(feature_name, 'spatial maps'))
        %% Spatial maps
        regressorsUnique = unique(regressors);%only 1 regressor per unique regressor
        out_files = write_composite_sig_effects(mancovanInfo.HInfo, regressorsUnique, result_files, ...
            prefix, outputDir, p_threshold, threshdesc);
        
        for nT = 1:length(out_files)
            tmpH = icatb_orth_views(out_files{nT}, 'structfile', structFile, 'image_values', image_values, 'convert_to_zscores', 'no', 'threshold', eps, ...
                'set_to_max_voxel', 1, 'fig_title', ['F-test ', regressorsUnique{nT}]);
            cH = findobj(tmpH, 'tag', 'compviewerColorbar');
            xlabel('-log_1_0(p-value)', 'parent', cH, 'Interpreter', 'tex')
            figH(end + 1).H = tmpH;
        end
        
    elseif (strcmpi(feature_name, 'timecourses spectra'))
        %% Timecourses spectra
        
        regressorsUnique = unique(regressors); %only 1 regressor per unique regressor
        for term = 1:length(regressorsUnique)
            figTitle = ['Univariate Results ', feature_name, ' F-test'];
            gH = icatb_getGraphics(figTitle, 'graphics', 'univariate_results_sepctra', 'on');
            
            load(fullfile(mancovanInfo.outputDir, result_files{1}), 'freq');
            msgStr = [feature_name, ' F-test ', regressorsUnique{term}, ' (p < ', num2str(mancovanInfo.display.p_threshold), ')'];
            
            [S, status] = gather_univariate_stats(mancovanInfo, nFeature, term, 0);
            if (~status)
                delete(gH);
                disp(['Feature ',feature_name, ': No ', msgStr]);
                continue;
            end
            
            sh = subplot(1, 1, 1);
            axesH = sh;
            cmap1 = coldhot;
            S.c = msgStr;
            vars1 = {S, freq, mancovanInfo.comps, mancovanInfo.display.p_threshold, mancovanInfo.display.threshdesc, cmap1, 'Component', ...
                'Frequency (Hz)', sh};
            
            ch1 = plotUnivStats (vars1{:});
            figH(end + 1).H = gH;
            
        end
        
    elseif (strcmpi(feature_name, 'fnc correlations'))
        %% FNC correlations
        regressorsUnique = unique(regressors); %only 1 regressor per unique regressor
        for term = 1:length(regressorsUnique)
            figTitle = ['Univariate Results ', feature_name];
            [S, status] = gather_univariate_stats(mancovanInfo, nFeature, term, 0, 1);
            msgStr = [feature_name, ' F-test ', regressorsUnique{term}, ' (p < ', num2str(mancovanInfo.display.p_threshold), ')'];
            if (~status)
                delete(gH);
                disp(['Feature ',feature_name, ': No ', msgStr]);
                continue;
            end
            
            gH = icatb_getGraphics(figTitle, 'graphics', 'univariate_results_fnc', 'on');
            
            nRows = 1;
            
            sh = subplot(nRows, 1, 1);
            axesH = sh;
            cmap1 = coldhot;
            S.c = msgStr;
            
            M = S.logp;
            M = icatb_vec2mat(M, 1);
            S2 = S;
            S2.p = S.p;
            S2.logp = M;
            vars = {S2, mancovanInfo.comps, mancovanInfo.comps, mancovanInfo.display.p_threshold, mancovanInfo.display.threshdesc, ...
                cmap1, 'Component', 'Component', sh, 1};
            ch1 = plotUnivStats (vars{:});
            %FNC_MAT = S.logp;
            
            clear vars;
            axis(sh, 'square');
            
            figH(length(figH) + 1).H = gH;
            
            if (display_connectogram)
                try
                    % Show connectogram
                    gH = icatb_plot_connectogram([], comp_network_names, 'C', M, 'image_file_names', mancovanInfo.userInput.compFiles, ...
                        'threshold', mancovanInfo.display.t_threshold, 'convert_to_zscores', 'yes', 'colorbar_label', ...
                        '-log10(p)', 'title', msgStr);
                    figH(length(figH) + 1).H = gH;
                catch
                end
            end
            
            
            if (length(result_files) > 1)
                
                
                load(fullfile(mancovanInfo.outputDir, result_files{2}), 'low_inds');
                
                % network average
                [S, status] = gather_univariate_stats(mancovanInfo, nFeature, term, 0, 2, comp_network_names(:,1));
                %[S, status] = gather_univariate_stats(mancovanInfo, nF, covariatesToPlot{nCov}, timeNo, 2, comp_network_names(:,1));
                msgStr = [feature_name, ' F-test ' , regressorsUnique{term}, ' (p < ', num2str(mancovanInfo.display.p_threshold), ') network averaged'];
                if (~status)
                    %delete(gH);
                    disp(['Feature ',feature_name, ': No ', msgStr]);
                    continue;
                end
                
                gH = icatb_getGraphics(figTitle, 'graphics', 'univariate_results_fnc', 'on');
                
                nRows = 1;
                
                
                sh = subplot(nRows, 1, 1);
                axesH = sh;
                cmap1 = coldhot;
                S.c = msgStr;
                M = S.logp;
                M = squeeze(icatb_low2Mat(M, length(comp_network_names(:,1)), low_inds));
                S2 = S;
                S2.p = S.p;
                S2.logp = M;
                vars = {S2, 1:length(comp_network_names(:,1)), comp_network_names(:,1), mancovanInfo.display.p_threshold, ...
                    mancovanInfo.display.threshdesc, cmap1, 'Component', 'Component', sh, 1};
                ch1 = plotUnivStats (vars{:});
                
                clear vars;
                axis(sh, 'square');
                
                figH(length(figH) + 1).H = gH;
                
                
                if (display_connectogram)
                    try
                        % Show connectogram
                        gH = icatb_plot_connectogram([], comp_network_names,'C', M, 'image_file_names', mancovanInfo.userInput.compFiles, ...
                            'iscomposite', 1, 'exclude_zeros', 0, ...
                            'comp_labels', comp_network_names(:, 1), 'threshold', mancovanInfo.display.t_threshold, ...
                            'convert_to_zscores', 'yes', 'title', msgStr, 'colorbar_label', '-log10(p)', ...
                            'cmap', coldhot);
                        figH(length(figH) + 1).H = gH;
                    catch
                    end
                end
            end
            % end for domain average
            
        end
        % for terms
        
        
    else
        %% FNC lag
        figTitle = ['Univariate Results ', feature_name];
        regressorsUnique = unique(regressors); %only 1 regressor per unique regressor
        for term = 1:length(regressorsUnique)
            [S, status] = gather_univariate_stats(mancovanInfo, nFeature, term, 0, 1:2, '', 1);
            
            %[S, status] = gather_univariate_stats(mancovanInfo, nF, covariatesToPlot{nCov}, timeNo, 1:2, '', 1);
            
            disp('p-values not surpassing the threshold are excluded from both lag and correlations ...');
            
            %chk_good_inds = (prod(S.p ~= 0));
            chk_good_inds = (isnan(S.p(1,:)) | isnan(S.p(2,:))) == 0;
            if (isempty(find(chk_good_inds == 1)))
                status = 0;
            end
            
            S.p(:, chk_good_inds == 0) = 0;
            S.logp(:, chk_good_inds == 0) = 0;
            
            msgStr = [feature_name, ' F-test ', regressorsUnique{term}, ' (p < ', num2str(mancovanInfo.display.p_threshold), ')'];
            
            if (~status)
                %delete(gH);
                disp(['Feature ',feature_name, ': No ', msgStr]);
                continue;
            end
            
            
            for nFncPlots = 1:size(S.logp, 1)
                gH = icatb_getGraphics(figTitle, 'graphics', 'univariate_results_fnc', 'on');
                load icatb_colors coldhot;
                %coldhot = coldhot_sensitive;
                nRows = 1;
                %for nCov = 1:length(covariatesToPlot)
                
                plot_arrows = 0;
                colorbarTitle = '-log10(p) (corr)';
                if (nFncPlots == 2)
                    colorbarTitle = '-log10(p) (Lag)';
                    msgStr = [feature_name, ' F-test ', regressorsUnique{term}, ' (p < ', num2str(mancovanInfo.display.p_threshold), ')'];
                    plot_arrows = 1;
                end
                
                
                sh = subplot(nRows, 1, 1);
                axesH = sh;
                cmap1 = coldhot;
                S.c = msgStr;
                M = S.logp(nFncPlots, :) ;
                M = icatb_vec2mat(M, 1);
                S2 = S;
                S2.p = S.p(nFncPlots, :);
                S2.logp = M;
                vars = {S2, mancovanInfo.comps, mancovanInfo.comps, mancovanInfo.display.p_threshold, mancovanInfo.display.threshdesc, ...
                    cmap1, 'Component', 'Component', sh, 1};
                ch1 = plotUnivStats (vars{:});
                
                clear vars;
                axis(sh, 'square');
                
                
                figH(length(figH) + 1).H = gH;
                
                if (display_connectogram)
                    try
                        % Show connectogram
                        gH = icatb_plot_connectogram([], comp_network_names, 'C', M, 'image_file_names', mancovanInfo.userInput.compFiles, ...
                            'threshold', mancovanInfo.display.t_threshold, ...
                            'convert_to_zscores', 'yes', 'colorbar_label', colorbarTitle, 'title', msgStr);
                        figH(length(figH) + 1).H = gH;
                    catch
                    end
                end
                
            end
            
        end
        
    end
end
% end for features

function outfiles = write_composite_sig_effects(V, terms, result_files, prefix, outputDir, thresh, threshdesc)
%% Write composite files for covariates of interest
%

outfiles = repmat({''}, 1, length(terms));

for nT = 1:length(terms)
    
    term = terms{nT};
    
    B = zeros(V(1).dim(1:3));
    
    tmpPath = fileparts(result_files{1});
    
    for ii = 1:length(result_files)
        
        load(fullfile(outputDir, result_files{ii}), 'UNI', 'mask_ind');
        tmp_mask_ind = mask_ind;
        
        TEMP = zeros(size(B));
        
        ps = UNI.FstatInfo.pANCOVAN(nT, :);
        fs = UNI.FstatInfo.FANCOVAN(nT, :);
        [p_masked, ps]  = get_sig_pvalues(ps, thresh, threshdesc);
        bad_inds = ps > p_masked;
        tmp_mask_ind(bad_inds) = [];
        ps = ps(bad_inds == 0);
        ps = -log10(ps + eps);
        TEMP(tmp_mask_ind) = ps;
        
        gmask = find(abs(TEMP) > abs(B));
        B(gmask) = TEMP(gmask);
        clear I mIND TEMP
    end
    
    
    V.fname = fullfile(outputDir, tmpPath, [prefix, '_sm_', term, '_fstat.nii']);
    
    V.n(1) = 1;
    icatb_write_vol(V, B);
    outfiles{nT} = V.fname;
    
end


function C = plotUnivStats(T, f, comps, thresh, threshdesc, cmap, yLabelStr, xLabelStr, axesH, fnc)

icatb_defaults;
global FONT_COLOR;

A = T.logp;
ncomp = size(A,1);
A(isnan(A)) = 0;

pv = T.p(find(isnan(T.p) == 0));
% if (strcmpi(threshdesc, 'fdr'))
%     p_masked  = icatb_fdr( pv, thresh);
% else
%     p_masked = thresh;
% end

[p_masked, pv]  = get_sig_pvalues(pv, thresh, threshdesc);

if (strcmpi(threshdesc, 'bhfdr'))
    %disp(p_masked)
    fdrlim = -log10(p_masked);
end

if (exist('fnc', 'var') && fnc)
    M = imagesc(1:length(comps),1:length(comps),A);
    set(axesH,'XTick',1:length(comps),'XTickLabel',T.y, 'TickDir', 'out')
    
else
    M = imagesc(f,1:length(comps),A);
end

CLIM = [-max(abs(A(:))) max(abs(A(:)))];

if (~isempty(find(CLIM ~= 0)))
    set(axesH, 'CLIM', CLIM);
end

set(axesH,'YTick',1:length(comps),'YTickLabel',T.y, 'TickDir', 'out')
%set(gca, 'CLim', [log10(eps) -log10(eps)])

colormap(cmap);
C = colorbar;
yt = get(C, 'YTick');
if (strcmpi(threshdesc, 'bhfdr'))
    set(C, 'YTick', sort([yt, -fdrlim fdrlim]));
end
ylabel(C, '-log_1_0(p-value)', 'Interpreter', 'tex');
if (length(find(A==0)) == numel(A))
    T.c = ['No ', T.c];
end
title(T.c);
ylabel(yLabelStr, 'parent', axesH);
xlabel(xLabelStr, 'parent', axesH);

set(axesH, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);

if (icatb_get_matlab_version >= 2014)
    set(C, 'color', FONT_COLOR);
end


function [p_masked, p]  = get_sig_pvalues(p, thresh, criteria)
% apply fdr correction

p_masked = thresh;
if (strcmpi(criteria, 'mafdr'))
    p = mafdr(p);
elseif(strcmpi(criteria, 'bhfdr'))
    p_masked = icatb_fdr(p, thresh);
end


function [S, status] = gather_univariate_stats(mancovanInfo, index, term, timeNo, resultsNum, LabelStr, fillNaN)
%% Gather univariate stats
%

result_files = mancovanInfo.outputFiles(index).filesInfo.result_files;
outputDir = mancovanInfo.outputDir;
thresh = mancovanInfo.display.p_threshold;
featureName = mancovanInfo.outputFiles(index).feature_name;
if (~exist('resultsNum', 'var'))
    resultInds = 1:length(result_files);
else
    resultInds = resultsNum;
end

if (~exist('fillNaN', 'var'))
    fillNaN = 0;
end

countR = 0;
for ii = resultInds
    countR = countR + 1;
    if (timeNo > 0)
        load(fullfile(outputDir, result_files{ii}), 'time', 'comp_number');
        UNI = time.UNI{timeNo};
    else
        load(fullfile(outputDir, result_files{ii}), 'UNI', 'comp_number');
    end
    regressors = UNI.FstatInfo.regressors;
    %mIND = strmatch(term, UNI.tests, 'exact');
    if (strcmpi(featureName, 'timecourses spectra'))
        load(fullfile(outputDir, result_files{ii}), 'freq');
        xdim = length(freq);
    elseif  (strcmpi(featureName, 'fnc correlations'))
        tmp_fnc_name = fullfile(outputDir, result_files{ii});
        varsInFNCFile = whos('-file', tmp_fnc_name);
        if (~isempty(strmatch('fnc_corrs', cellstr(char(varsInFNCFile.name)), 'exact')))
            load(tmp_fnc_name, 'fnc_corrs');
        end
        
        if (exist('fnc_corrs', 'var'))
            xdim = size(fnc_corrs, 2);
        else
            load(tmp_fnc_name, 'UNI');
            xdim = size(UNI.t{1}, 2);
        end
    else
        tmp_fnc_name = fullfile(outputDir, result_files{ii});
        varsInFNCFile = whos('-file', tmp_fnc_name);
        if (~isempty(strmatch('fnc_values', cellstr(char(varsInFNCFile.name)), 'exact')))
            load(tmp_fnc_name, 'fnc_values');
        end
        
        if (exist('fnc_values', 'var'))
            xdim = size(fnc_values, 2);
        else
            load(tmp_fnc_name, 'UNI');
            xdim = size(UNI.t{1}, 2);
        end
        
        %load(fullfile(outputDir, result_files{ii}), 'fnc_corrs');
        %xdim = size(fnc_corrs, 2);
        %clear fnc_corrs;
    end
    
    if (countR == 1)
        if (fillNaN)
            STATS = NaN(length(resultInds), xdim);
        else
            STATS = zeros(length(resultInds), xdim);
        end
        STATS2 = STATS;
    end
    
    ps = UNI.FstatInfo.pANCOVAN(term, :);
    fs = UNI.FstatInfo.FANCOVAN(term, :);
    [p_masked, ps]  = get_sig_pvalues(ps, thresh, mancovanInfo.display.threshdesc);
    %         if (strcmpi(mancovanInfo.display.threshdesc, 'fdr'))
    %             p_masked = icatb_fdr(ps, thresh);
    %         else
    %             p_masked = thresh;
    %         end
    good_inds = ps < p_masked;
    STATS(countR, good_inds) = -log10(ps(good_inds) + eps);
    STATS2(countR, good_inds) = ps(good_inds);
    conname = regressors{term};
    
end


if (~exist('conname', 'var'))
    conname = term;
end


status = 1;
if (length(find(STATS==0)) == numel(STATS))
    status = 0;
end

S.logp = STATS;
S.p = STATS2;
S.x = term;
S.c = conname;

if (exist('LabelStr', 'var') && ~isempty(LabelStr))
    S.y = LabelStr;
else
    S.y = cellstr(num2str(mancovanInfo.comps(:)))';
end


function fname = capitalize(fname)
% Convert first letter to upper case

fname = [upper(fname(1)), fname(2:end)];