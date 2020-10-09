function icatb_display_spatial_chronnectome(schronn_param_file)
%% Display spatial chronnectome results
%

icatb_defaults;

if (~exist('schronn_param_file','var'))
    schronn_param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Spatial Chronnectome Parameter File', 'filter', '*schronn.mat');
    drawnow;
end

if (isempty(schronn_param_file))
    error('Spatial chronnectome parameter file not selected');
end

load(schronn_param_file);

outputDir = fileparts(schronn_param_file);
if (isempty(outputDir))
    outputDir = pwd;
end


%sdhInfo.outputDir = outputDir;

results_files = schronnInfo.postprocess.results_files;

% For each functional domain display centroids
countIm = 0;
graphicsH = [];
for nF = 1:length(results_files)
    output_file = results_files{nF};
    load(fullfile(outputDir, output_file));
    for ii = 1:length(clusterInfo.clusterFiles)
        countIm = countIm + 1;
        tmpStr = sprintf('%d (%d%%)', sum(clusterInfo.IDXall == ii),  round(100*sum(clusterInfo.IDXall == ii)/length(clusterInfo.IDXall)));
        titleStr = ['Comp ',  icatb_returnFileIndex(schronnInfo.comps(nF)), ' Cluster Centroid (', num2str(ii), ')'];
        titleStr= [titleStr, ' ', tmpStr];
        if (countIm == 1)
            disp_parameters = icatb_image_viewer(fullfile(outputDir, clusterInfo.clusterFiles{ii}), 'usegui', 1, 'labels', {titleStr});
        else
            icatb_image_viewer(fullfile(outputDir, clusterInfo.clusterFiles{ii}), 'usegui', 0, 'labels', {titleStr}, ...
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
