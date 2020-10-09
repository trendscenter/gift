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
            
        else
            fig_title = ['Cluster Centroid (', num2str(ii), ')'];
            H =  icatb_getGraphics(fig_title, 'graphics', 'ROI-ROI Cluster States', 'on');
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
        end
        
        drawnow;
        graphicsH(end+1).H = gcf;
    end
    % end for centroids
end
% end for rois

drawnow;
icatb_plotNextPreviousExitButtons(graphicsH);

