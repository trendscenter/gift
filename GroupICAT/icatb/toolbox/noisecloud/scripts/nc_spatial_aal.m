function feature = nc_spatial_aal(spatialMap, threshold);
    global atlas;
        % AAL ATLAS FEATURES;
    % Requires aal2mni image registered to spatialMap space, see;
    % INSTRUCTIONS.txt for details;
        % Labels are included in AALlabels.txt and AAL_labels.mat, but will;
    % be read from database for feature extraction:;
    %   Precentral_L,Precentral_R,Frontal_Sup_L,Frontal_Sup_R,Frontal_Sup_Orb_L,Frontal_Sup_Orb_R,Frontal_Mid_L,Frontal_Mid_R,Frontal_Mid_Orb_L,Frontal_Mid_Orb_R,Frontal_Inf_Oper_L,Frontal_Inf_Oper_R,Frontal_Inf_Tri_L,Frontal_Inf_Tri_R,Frontal_Inf_Orb_L,Frontal_Inf_Orb_R,Rolandic_Oper_L,Rolandic_Oper_R,Supp_Motor_Area_L,Supp_Motor_Area_R,Olfactory_L,Olfactory_R,Frontal_Sup_Medial_L,Frontal_Sup_Medial_R,Frontal_Med_Orb_L,Frontal_Med_Orb_R,Rectus_L,Rectus_R,Insula_L,Insula_R,Cingulum_Ant_L,Cingulum_Ant_R,Cingulum_Mid_L,Cingulum_Mid_R,Cingulum_Post_L,Cingulum_Post_R,Hippocampus_L,Hippocampus_R,ParaHippocampal_L,ParaHippocampal_R,Amygdala_L,Amygdala_R,Calcarine_L,Calcarine_R,Cuneus_L,Cuneus_R,Lingual_L,Lingual_R,Occipital_Sup_L,Occipital_Sup_R,Occipital_Mid_L,Occipital_Mid_R,Occipital_Inf_L,Occipital_Inf_R,Fusiform_L,Fusiform_R,Postcentral_L,Postcentral_R,Parietal_Sup_L,Parietal_Sup_R,Parietal_Inf_L,Parietal_Inf_R,SupraMarginal_L,SupraMarginal_R,Angular_L,Angular_R,Precuneus_L,Precuneus_R,Paracentral_Lobule_L,Paracentral_Lobule_R,Caudate_L,Caudate_R,Putamen_L,Putamen_R,Pallidum_L,Pallidum_R,Thalamus_L,Thalamus_R,Heschl_L,Heschl_R,Temporal_Sup_L,Temporal_Sup_R,Temporal_Pole_Sup_L,Temporal_Pole_Sup_R,Temporal_Mid_L,Temporal_Mid_R,Temporal_Pole_Mid_L,Temporal_Pole_Mid_R,Temporal_Inf_L,Temporal_Inf_R,Cerebelum_Crus1_L,Cerebelum_Crus1_R,Cerebelum_Crus2_L,Cerebelum_Crus2_R,Cerebelum_3_L,Cerebelum_3_R,Cerebelum_4_5_L,Cerebelum_4_5_R,Cerebelum_6_L,Cerebelum_6_R,Cerebelum_7b_L,Cerebelum_7b_R,Cerebelum_8_L,Cerebelum_8_R,Cerebelum_9_L,Cerebelum_9_R,Cerebelum_10_L,Cerebelum_10_R,Vermis_1_2,Vermis_3,Vermis_4_5,Vermis_6,Vermis_7,Vermis_8,Vermis_9,Vermis_10;
    % Get all indices in atlas, each is a different region;
    atlas_index = unique(atlas(:));
    % Get rid of atlas_index(1), which is the index for 0;
    atlas_index(1) = [];
    feature = zeros(size(atlas_index,1),1);

    % salman 20200525 
    % without thresholding every SM returns the same AAL features
    if ~isnan(threshold)
        spatialMap( abs( spatialMap(:) ) < threshold ) = 0;
    end

        % For each region in the AAL mask (with a unique ID), create a mask;
    % for the region, and multiply by the current image to get voxels;
    % with activation.  Sum this image to get a count, and save to vector;
    for j=1:length(atlas_index);
        feature(j) = length(find(spatialMap(atlas == atlas_index(j)~= 0)));
    end;
end
