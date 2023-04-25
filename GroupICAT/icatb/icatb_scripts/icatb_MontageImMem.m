function [im3] = icatb_MontageImMem(DIstruImComponent, struImBgTemplate, bConvertToZ, dThreshold, drowSlices_in_mm)
% plots im from mem
    % struImComponent - single component loaded with icatb_spm_vol (often icasig component)
    % struImBgTemplate - structural volume
    % bConvertToZ: 0 or 1 depending on im scaled to z score
    % dThreshold: level of theshold
    % drowSlices_in_mm: pos of slice cuts. e.g. [-40   -30   -20   -10     0    10    20    30    40    50    60    70]
    % returning im3 which is the image
    %example:
    % im3 = funMontageImMem(stru_icasig,struImBgTemplate, 1, 1, [-40 -30 -20 -10 0 10 20 30 40 50 60 70]);

    % Initiation
    anatomical_view = 'axial';    

    % Set color map  
    load icatb_colors coldhot cold hot;
    cmap = coldhot(1:4:end, :); % positive and negative
    cmap = [cmap; gray(size(cmap, 1))];

    % Axes
    gH = icatb_getGraphics('My Title', 'graphics', 'Image Viewer', 'on');
    axesH = axes('parent', gH, 'units', 'normalized', 'position', [0.1, 0.1, 0.8, 0.8]);
    
    %Resize to structural
    [images, coords, HInfo, slices_in_mm, text_left_right] = icatb_MontageImMemResize(struImBgTemplate, DIstruImComponent, anatomical_view,  drowSlices_in_mm);

    % Get Image
    irowImNumAndSize = size(images);
    DIM=irowImNumAndSize(2:4);

    struImBgTemplate1d = reshape(images(1,:,:,:), 1, prod(DIM));
    
    returnValue = 1;
    im3MyComponentMod=reshape(images(2,:,:,:), 1, prod(DIM));
    im3MyComponentModZ = icatb_applyDispParameters(im3MyComponentMod, bConvertToZ, returnValue, dThreshold, [1, prod(DIM)], HInfo);
    
    [DIstruImComponentZ,  maxICAIM, minICAIM, minInterval, maxInterval] = icatb_check_overlayComp(im3MyComponentModZ, struImBgTemplate1d, 1, 1);
    DIstruImComponentZ = reshape(DIstruImComponentZ, DIM);
    [im3, numImagesX, numImagesY, drowSlices_in_mm_new] = icatb_returnMontage(DIstruImComponentZ, [], DIM, HInfo.VOX, drowSlices_in_mm);
    colormap(cmap);
    ImageAxis = image(im3, 'parent', axesH, 'CDataMapping', 'scaled');
    set(axesH, 'clim',  [minInterval, 2*maxInterval]); % set the axis positions to the specified
    axis(axesH, 'off');
    axis(axesH, 'image')
    
end

