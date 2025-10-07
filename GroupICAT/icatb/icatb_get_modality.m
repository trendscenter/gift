function [modalityType, dataTitle, compSetFields] = icatb_get_modality
% Return modality type and data title

appName = 'group_ica_modality';

if isappdata(0, appName)
    group_ica_modality = getappdata(0, appName);
else
    group_ica_modality = 'fmri';
end

if isempty(group_ica_modality)
    group_ica_modality = 'fmri';
end

if strcmpi(group_ica_modality, 'fmri')
    modalityType = 'fMRI';
    dataTitle = 'Functional';
    compSetFields = {'ic', 'tc'};
elseif strcmpi(group_ica_modality, 'smri')
    modalityType = 'sMRI';
    dataTitle = 'Structural';
    compSetFields = {'ic', 'tc'};
    
elseif strcmpi(group_ica_modality, 'fnc')
    modalityType = 'FNC';
    dataTitle = 'FNC';
    compSetFields = {'ic', 'tc'};
elseif strcmpi(group_ica_modality, 'conn')
    modalityType = 'CONN';
    dataTitle = 'CONN';
    compSetFields = {'ic', 'none'};
else
    modalityType = 'EEG';
    dataTitle = 'EEG';
    compSetFields = {'timecourse', 'topography'};
end

