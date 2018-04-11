function graphicsH = icatb_plot_mancova_features(mancovanInfo)

icatb_defaults;
global UI_FONTNAME;
global UI_FS;

features = mancovanInfo.features;
comps = mancovanInfo.comps;
sm_inds = strmatch('spatial maps', lower(features), 'exact');
spectra_inds = strmatch('timecourses spectra', lower(features), 'exact');
fnc_inds = strmatch('fnc correlations', lower(features), 'exact');
fnc_inds_lag = strmatch('fnc correlations (lag)', lower(features), 'exact');
fnc_inds = [fnc_inds, fnc_inds_lag];

try
    freq_limits = mancovanInfo.freq_limits;
catch
    freq_limits = [0.1, 0.15];
end

try
    structFile = mancovanInfo.structFile;
catch
    structFile = fullfile(fileparts(which('gift.m')), 'icatb_templates', 'ch2bet.nii');
end

try
    t_threshold = mancovanInfo.t_threshold;
catch
    t_threshold = 1;
end

try
    imStr = mancovanInfo.image_values;
catch
    imStr = 'Positive';
end


plotSM = 0;
plotTC = 0;
if (~isempty(sm_inds))
    plotSM = 1;
    %sm_results = results{sm_inds};
    %tmap_files = cellstr(char(sm_results.tmap_file));
end

if (~isempty(spectra_inds))
    plotTC = 1;
    %spectra_results = results{spectra_inds};
end

load icatb_colors coldhot;
graphicsH = [];

if (plotSM || plotTC)
    %graphicsH = repmat(struct('H', []), 1, length(comps));
    %% Display ortho slices and spectra
    countF = 0;
    for nFiles = 1:length(comps)
        countF = countF + 1;
        
        if (plotTC)
            load(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(spectra_inds).filesInfo.result_files{nFiles}), 'spectra_tc', 'freq');
            if (~exist('spectra_tc', 'var') || isempty(spectra_tc))
                error('Please run setup features in order to view them.');
            end
            tc.data = spectra_tc;
            tc.xAxis = freq;
            tc.isSpectra = 1;
            tc.xlabelStr = 'Frequency (Hz)';
            tc.ylabelStr = 'Power';
            dynamicrange = zeros(1, size(tc.data, 1));
            fALFF = dynamicrange;
            for nS = 1:size(tc.data, 1)
                [dynamicrange(nS), fALFF(nS)] = icatb_get_spec_stats(tc.data(nS, :), tc.xAxis, freq_limits);
            end
            dynamicrange = mean(dynamicrange);
            fALFF = mean(fALFF);
            tc.titleStr = sprintf('Dynamic range: %0.3f, Power_L_F/Power_H_F: %0.3f', mean(dynamicrange), mean(fALFF));
        end
        
        if (plotTC && plotSM)
            H = icatb_orth_views(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(sm_inds).filesInfo.tmap_files{nFiles}), 'structfile', structFile, 'image_values', imStr, 'convert_to_zscores', 'no', 'set_to_max_voxel', 1, 'tc', tc, 'labels', ...
                ['T-map ', icatb_returnFileIndex(comps(nFiles)), ' (T >= ', num2str(t_threshold), ')'], 'fig_title', ['Features (Tmap ', icatb_returnFileIndex(comps(nFiles)), ' and Spectra)'], 'threshold', t_threshold);
        else
            if (plotSM)
                H = icatb_orth_views(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(sm_inds).filesInfo.tmap_files{nFiles}), 'structfile', structFile, 'image_values', imStr, 'convert_to_zscores', 'no', 'set_to_max_voxel', 1, 'labels', ...
                    ['T-map ', icatb_returnFileIndex(comps(nFiles)), ' (T >= ', num2str(t_threshold), ')'], 'fig_title', ['Features (Tmap ', icatb_returnFileIndex(comps(nFiles)), ')'], 'threshold', t_threshold);
            else
                tc.titleStr = sprintf('Comp %s Dynamic range: %0.3f, Power_L_F/Power_H_F: %0.3f', icatb_returnFileIndex(comps(nFiles)), mean(dynamicrange), mean(fALFF));
                tc.fig_title = ['Power Spectra ', icatb_returnFileIndex(comps(nFiles))];
                H = icatb_plot_spectra(tc);
                axisH = axes('Parent', H, 'position', [0 0 1 1], 'visible', 'off');
                xPos = 0.5; yPos = 0.97;
                titleColor = 'c';
                text(xPos, yPos, tc.fig_title, 'color', titleColor, 'fontweight', 'bold', 'fontsize', UI_FS + 2, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME, 'parent', axisH);
            end
        end
        
        graphicsH(length(graphicsH) + 1).H = H;
    end
    
end

clear spectra_tc

if (~isempty(fnc_inds))
    
    network_values = zeros(1, length(mancovanInfo.userInput.comp));
    for nV = 1:length(network_values)
        network_values(nV) = length(mancovanInfo.userInput.comp(nV).value);
    end
    network_names =  cellstr(char(mancovanInfo.userInput.comp.name));
    
    mycmap = jet(64);
    
    %varsIn = {'fnc_corrs', 'fnc_values'};
    
    for nFnc = 1:length(fnc_inds)
        
        if (~isempty(icatb_findstr(features{fnc_inds(nFnc)}, 'lag')))
            varsIn = 'fnc_values';
        else
            varsIn = 'fnc_corrs';
        end
        
        tmpFnc = load(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(fnc_inds(nFnc)).filesInfo.result_files{1}), varsIn);
        
        fnc_corrs = tmpFnc.(varsIn);
        
        compValues = mancovanInfo.comps(:);
        compValues = cellstr(num2str(compValues));
        
        M = icatb_vec2mat(icatb_z_to_r(squeeze(mean(fnc_corrs))), 1);
        CLIM = max(abs(M(:)));
        fig_title = ['Features ', features{fnc_inds(nFnc)}];
        axesTitle = 'FNC Correlations (Averaged over subjects)';
        isLag = ~isempty(icatb_findstr(features{fnc_inds(nFnc)}, 'lag'));
        if (isLag)
            axesTitle = 'FNC Correlations with lag (Averaged over subjects)';
        end
        
        
        gH = plotFNCResults(M, compValues, CLIM, mycmap, fig_title, axesTitle, network_names, network_values);
        
        
        %         if (length(network_names) == 1)
        %
        %             gH = icatb_plot_matrix(M, cellstr(num2str(mancovanInfo.comps(:))), cellstr(num2str(mancovanInfo.comps(:))), 'fig_title', fig_title, 'tag', 'FNC Correlations', 'title', ...
        %                 axesTitle, 'cmap', mycmap, 'clim', [-CLIM, CLIM], 'ylabel', 'Components', 'xlabel', 'Components');
        %
        %         else
        %
        %             gH = icatb_getGraphics(fig_title, 'graphics',  'FNC Correlations', 'on');
        %             set(gH, 'resize', 'on');
        %             axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
        %             icatb_plot_FNC(M, [-CLIM, CLIM], cellstr(num2str(mancovanInfo.comps(:))), (1:length(mancovanInfo.comps)), gH, axesTitle, axesH, ...
        %                 network_values, network_names);
        %             colormap(mycmap);
        %
        %         end
        
        %         axisH = axes('Parent', gH, 'position', [0 0 1 1], 'visible', 'off');
        %         xPos = 0.5; yPos = 0.97;
        %         titleColor = 'c';
        %         text(xPos, yPos, fig_title, 'color', titleColor, 'fontweight', 'bold', 'fontsize', UI_FS + 2, 'HorizontalAlignment', 'center', 'FontName', UI_FONTNAME, 'parent', axisH);
        
        graphicsH(length(graphicsH) + 1).H = gH;
        
        
        
        if (isLag)
            axesTitle = 'Lag in seconds (Averaged over subjects)';
            fig_title = axesTitle;
            tmpFnc = load(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(fnc_inds(nFnc)).filesInfo.result_files{2}), varsIn);
            fnc_corrs = tmpFnc.(varsIn);
            M = icatb_vec2mat((squeeze(mean(fnc_corrs))), 1);
            CLIM = max(abs(M(:)));
            gH = plotFNCResults(M, compValues, CLIM, mycmap, fig_title, axesTitle, network_names, network_values);
            graphicsH(length(graphicsH) + 1).H = gH;
        end
        
        
        
        if (~isLag)
            % Domain averaged fnc
            
            if (length(mancovanInfo.outputFiles(fnc_inds(nFnc)).filesInfo.result_files) == 2)
                
                axesTitle = 'FNC Correlations domain averaged (Averaged over subjects)';
                
                tmpFnc = load(fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(fnc_inds(nFnc)).filesInfo.result_files{2}), varsIn, 'low_inds');
                
                fnc_corrs = tmpFnc.(varsIn);
                lower_inds = tmpFnc.('low_inds');
                
                compValues = (1:length(network_names))';
                compValues = cellstr(num2str(compValues));
                
                M = squeeze(icatb_low2Mat(icatb_z_to_r(squeeze(mean(fnc_corrs))), length(compValues), lower_inds));
                
                %M = vec2matN(icatb_z_to_r(squeeze(mean(fnc_corrs))), lower_inds);
                CLIM = max(abs(M(:)));
                fig_title = ['Features ', features{fnc_inds(nFnc)}, ' domain averaged'];
                
                %network_valuesT = ones(1, length(network_names));
                gH = icatb_plot_matrix(M, network_names, network_names, 'fig_title', fig_title, 'tag', 'FNC Correlations', 'title', ...
                    axesTitle, 'cmap', mycmap, 'clim', [-CLIM, CLIM], 'ylabel', 'Correlations', 'xlabel', ' ');
                %gH = plotFNCResults(M, network_names, CLIM, mycmap, fig_title, axesTitle, network_names, network_valuesT);
                
                graphicsH(length(graphicsH) + 1).H = gH;
                
            end
            
        end
        
        
        
    end
    
end

clear fnc_corrs;



function gH = plotFNCResults(M, compValues, CLIM, mycmap, fig_title, axesTitle, network_names, network_values)


if (length(network_names) == 1)
    
    gH = icatb_plot_matrix(M, compValues, compValues, 'fig_title', fig_title, 'tag', 'FNC Correlations', 'title', ...
        axesTitle, 'cmap', mycmap, 'clim', [-CLIM, CLIM], 'ylabel', 'Components', 'xlabel', 'Components');
    
else
    
    gH = icatb_getGraphics(fig_title, 'graphics',  'FNC Correlations', 'on');
    set(gH, 'resize', 'on');
    axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    icatb_plot_FNC(M, [-CLIM, CLIM], compValues, (1:length(compValues)), gH, axesTitle, axesH, ...
        network_values, network_names);
    colormap(mycmap);
    
end