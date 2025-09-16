function icatb_display_sdh
%% Display spatial dynamics hierarchy results
%

sdh_param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Spatial Dynamics Parameter File', 'filter', '*sdh*info*mat');
drawnow;
if (isempty(sdh_param_file))
    error('SDH parameter file not selected');
end

load(sdh_param_file);

outputDir = fileparts(sdh_param_file);
if (isempty(outputDir))
    outputDir = pwd;
end

%sdhInfo.outputDir = outputDir;

% For each functional domain display centroids
graphicsH = [];
countIm = 0;
for nF = 1:length(sdhInfo.postprocess.outputFiles)
    output_file = fullfile(outputDir, sdhInfo.postprocess.outputFiles{nF});
    load(output_file);
    clusterInfo = domainInfo.clusterInfo;
    for ii = 1:length(domainInfo.clusterFiles)
        countIm = countIm + 1;
        tmpStr = sprintf('%d (%d%%)', sum(clusterInfo.IDXall == ii),  round(100*sum(clusterInfo.IDXall == ii)/length(clusterInfo.IDXall)));
        titleStr = ['Domain ',  sdhInfo.comp_network_names{nF, 1}, ' Cluster Centroid (', num2str(ii), ')'];
        titleStr= [titleStr, ' ', tmpStr];
        if (countIm == 1)
            disp_parameters = icatb_image_viewer(fullfile(outputDir, domainInfo.clusterFiles{ii}), 'usegui', 1, 'labels', {titleStr});
        else
            icatb_image_viewer(fullfile(outputDir, domainInfo.clusterFiles{ii}), 'usegui', 0, 'labels', {titleStr}, ...
                'structfile', disp_parameters.structFile, 'image_values', disp_parameters.image_values, 'threshold', disp_parameters.threshold, ...
                'convert_to_zscores', disp_parameters.convert_to_z, 'slices_in_mm', disp_parameters.slices_in_mm, ...
                'slice_plane', disp_parameters.slice_plane, 'display_type', disp_parameters.display_type, 'isComposite', disp_parameters.isComposite);
        end
        drawnow;
        graphicsH(end+1).H = gcf;
    end
    
end

drawnow;
icatb_plotNextPreviousExitButtons(graphicsH);

