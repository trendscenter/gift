function icatb_display_dfc_roi
%% Display dFC ROI results
%

icatb_defaults;
global UI_FONTNAME;
global UI_FS;
global FONT_COLOR;

dfc_roi_param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select dFC ROI parameter file', 'filter', '*dfc*roi.mat');
drawnow;
if (isempty(dfc_roi_param_file))
    error('dFC ROI parameter file not selected');
end

load(dfc_roi_param_file);

outputDir = fileparts(dfc_roi_param_file);
if (isempty(outputDir))
    outputDir = pwd;
end

roi_labels_sel = dfcRoiInfo.postprocess.roi_labels_sel;
analysisType = dfcRoiInfo.userInput.analysisType;


% For each functional domain display centroids
graphicsH = [];
countIm = 0;
for nF = 1:length(dfcRoiInfo.postprocess.outputFiles)
    output_file = fullfile(outputDir, dfcRoiInfo.postprocess.outputFiles{nF});
    load(output_file);
    %clusterInfo = domainInfo.clusterInfo;
    for ii = 1:clusterInfo.num_clusters
        countIm = countIm + 1;
        tmpStr = sprintf('%d (%d%%)', sum(clusterInfo.IDXall == ii),  round(100*sum(clusterInfo.IDXall == ii)/length(clusterInfo.IDXall)));
        
        if (strcmpi(analysisType, 'roi-voxel'))
            
            titleStr = ['ROI ',  roi_labels_sel{nF}, ' Cluster Centroid (', num2str(ii), ')'];
            titleStr= [titleStr, ' ', tmpStr];
            if (countIm == 1)
                disp_parameters = icatb_image_viewer(fullfile(outputDir, clusterInfo.clusterFiles{ii}), 'usegui', 1, 'labels', {titleStr});
            else
                icatb_image_viewer(fullfile(outputDir, clusterInfo.clusterFiles{ii}), 'usegui', 0, 'labels', {titleStr}, ...
                    'structfile', disp_parameters.structFile, 'image_values', disp_parameters.image_values, 'threshold', disp_parameters.threshold, ...
                    'convert_to_zscores', disp_parameters.convert_to_z, 'slices_in_mm', disp_parameters.slices_in_mm, ...
                    'slice_plane', disp_parameters.slice_plane, 'display_type', disp_parameters.display_type, 'isComposite', disp_parameters.isComposite);
            end
            graphicsH(end+1).H = gcf;
            drawnow;
        else
            fig_title = ['Cluster Centroid (', num2str(ii), ')'];
            H =  icatb_getGraphics(fig_title, 'graphics', 'ROI-ROI Cluster States', 'on');
            figure(H);
            colormap(jet);
            sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
            CLIM = max(abs(clusterInfo.Call(:)));
            CLIM = [-CLIM, CLIM];
            tmp = icatb_vec2mat(clusterInfo.Call(ii, :), 1);
            titleStr = tmpStr;
            network_vals = ones(length(roi_labels_sel), 1);
            %cellstr(num2str((1:length(roi_labels_sel))'))
            icatb_plot_FNC(tmp, CLIM, roi_labels_sel, (1:length(roi_labels_sel)), H, 'Correlations (z)', ...
                sh(1), network_vals, roi_labels_sel);
            title(titleStr, 'parent', sh(1), 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            axis(sh, 'square');
            set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            graphicsH(end+1).H = H;
            drawnow;
        end
    end
    % end for centroids
end
% end for rois

if ~(strcmpi(analysisType, 'roi-voxel'))
    % Variables needed for summary statistics
    iSub = size(clusterInfo.states,1); % Number of subjects
    iSess = size(clusterInfo.states,2); % Number of sessions
    Nwin = length(clusterInfo.IDXall) / iSub;
    aIND = reshape(clusterInfo.IDXall, iSub, Nwin);
    aIND = aIND';
    
    % Generate the summary statistics
    matFractStateTime = zeros(iSub, clusterInfo.num_clusters);
    matTransitions = zeros(iSub, clusterInfo.num_clusters, clusterInfo.num_clusters);
    matMeanDwellTime = zeros(iSub, clusterInfo.num_clusters);
    matNumTransitions = zeros(iSub, 1);
    for ii = 1:iSub
        [FRii, TMii, MDTii, NTii] = icatb_dfnc_statevector_stats(aIND(:,ii), clusterInfo.num_clusters);
        matFractStateTime(ii,:) = FRii;
        matTransitions(ii,:,:) = TMii;
        matMeanDwellTime(ii,:) = MDTii;
        matNumTransitions(ii) = NTii;
    end
    
    % save the summary statistics
    README_icatb={'Score Summary Information'; ...
                'Type following commands to get statistacal reports'; ...
                'mean(matFractStateTime(:,:)) % returns the average fractal state time for the group'; ...
                'squeeze(mean(matTransitions(:,:,:))) % returns the mean of the transition matrix' ; ...
                'mean(matMeanDwellTime) % returns the mean of the subject mean dwell times' ; ...
                'matNumTransitions % returns the number of transitions per subject'};
    sDir = fileparts(dfc_roi_param_file);
    save(fullfile(sDir, [dfcRoiInfo.prefix '_dfc_display_sum_scores.mat']), 'README_icatb', 'matNumTransitions', 'matFractStateTime', 'matTransitions', 'matMeanDwellTime');
    disp(['dFC summary statistics saved to ' fullfile(sDir, [dfcRoiInfo.prefix '_dfc_display_sum_scores.mat'])]);
    
    % Plot the summary statistics
    figVisible = true;
    graphicsH(end+1).H = icatb_getGraphics('State Vector Stats', 'graphics', 'dfnc_summary4', figVisible);
    colormap(jet);
    
    icatb_dfnc_plot_statevector_stats(clusterInfo.num_clusters, matFractStateTime, matTransitions, matMeanDwellTime, matNumTransitions, graphicsH(end).H);
    set(findobj(graphicsH(end).H, 'type', 'axes'), 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
    
    try
        set(findobj(graphicsH(end).H, 'type', 'colorbar'), 'color', FONT_COLOR);
    catch
    end
    drawnow;
    
    % HTML for future use
    if (nargout == 1)
       titleStr = 'State Vector Stats';
       outFile = [dfcRoiInfo.prefix '_state_vector_stats.jpg'];
       printFile(graphicsH(end).H, fullfile(sDir, outFile));
       results(end + 1).file = outFile;
       results(end).text = 'Plot shows frequency of each cluster followed by mean dwell time in windows followed by mean of state transition matrix across subjects';
    end
    
    %% State vector
    if (~isfield(clusterInfo, 'states'))
        % subjects x sessions x windows
        states = reshape(clusterInfo.IDXall, iSess, iSub, Nwin);
        states = permute(states, [2, 1, 3]);
    else
        states = clusterInfo.states;
        states = reshape(states, iSub, iSess, Nwin);
    end
    
    numStateRows = ceil(sqrt(iSub));
    numStateCols = ceil(iSub/numStateRows);
    for nSess = 1:iSess
        graphicsH(end + 1).H = icatb_getGraphics(['State vector for session ', num2str(nSess)], 'graphics', ['dfnc_summary', num2str(4 + nSess)], figVisible);
        for nSub = 1:iSub
            sh = subplot(numStateRows, numStateCols, nSub);
            plot(squeeze(states(nSub, nSess, :)), 'm', 'parent', sh);
            title(['Sub ', num2str(nSub)], 'parent', sh, 'fontname', UI_FONTNAME, 'fontsize', 7);
            axis(sh, [0, Nwin, 0, clusterInfo.num_clusters+1]);
            if (nSub ~= 1)
                set(sh, 'XTick', []);
                set(sh, 'YTick', []);
            end
        end
        set(findobj(graphicsH(end).H, 'type', 'axes'), 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
        drawnow;
        if (nargout == 1) % HTML for future use
            titleStr = ['States for session ', num2str(nSess)];
            outFile = [dfcRoiInfo.prefix, '_state_vector_session', num2str(nSess), '.jpg'];
            printFile(graphicsH(end).H, fullfile(sDir, outFile));
            results(end + 1).file = outFile;
            results(end).text = ['State vector for session ', num2str(nSess)];
        end
    end
        
    % HTML for future use
    if (nargout == 1)
       for nH = 1:length(graphicsH)
           close(graphicsH(nH).H);
       end
       varargout{1} = results;
       return;
    end
end

drawnow;
icatb_plotNextPreviousExitButtons(graphicsH);
drawnow;


