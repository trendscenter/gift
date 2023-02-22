function default_mask_opts = icatb_default_mask_opts
%% Default mask options
%

icatb_defaults;
global DEFAULT_MASK_OPTION;
global REMOVE_CONSTANT_VOXELS;
global DEFAULT_MASK_SBM_MULTIPLIER;

modalityType = icatb_get_modality;

dm_sbm_mult = DEFAULT_MASK_SBM_MULTIPLIER;
if (isempty(dm_sbm_mult))
    dm_sbm_mult = 0.01;
end

removeConstVoxels = 0;
if (~isempty(REMOVE_CONSTANT_VOXELS))
    removeConstVoxels = REMOVE_CONSTANT_VOXELS;
end

defaultMaskOption = 'first_file';
if (~isempty(DEFAULT_MASK_OPTION))
    defaultMaskOption = DEFAULT_MASK_OPTION;
end

default_mask_opts = [];

if (strcmpi(modalityType, 'fmri'))
    default_mask_opts = struct('use_all_files', strcmpi(defaultMaskOption, 'all_files'), 'remove_constant_voxels', removeConstVoxels);
elseif (strcmpi(modalityType, 'smri'))
    default_mask_opts = struct('dm_mult', dm_sbm_mult);
end