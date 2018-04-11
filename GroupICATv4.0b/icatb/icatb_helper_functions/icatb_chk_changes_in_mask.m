function icatb_chk_changes_in_mask(sesInfo)
%% Check changes in default mask
%

modalityType = icatb_get_modality;

if (~isfield(sesInfo.userInput, 'mask_ind') || isempty(sesInfo.userInput.mask_ind))
    setappdata(0, 'create_mask_gica', 1);
    return;
end

if (isfield(sesInfo.userInput, 'default_mask_opts'))
    oldOpts = sesInfo.userInput.default_mask_opts;
    newOpts = icatb_default_mask_opts;
    if (isfield(sesInfo.userInput, 'maskFile') && isempty(sesInfo.userInput.maskFile))
        if (isfield(sesInfo.userInput, 'default_mask_opts'))
            if (strcmpi(modalityType, 'fmri'))
                if (oldOpts.use_all_files ~= newOpts.use_all_files)
                    setappdata(0, 'create_mask_gica', 1);
                    return;
                end
                if (oldOpts.remove_constant_voxels ~= newOpts.remove_constant_voxels)
                    setappdata(0, 'create_mask_gica', 1);
                    return;
                end
            elseif (strcmpi(modalityType, 'smri'))
                if (oldOpts.dm_mult ~= newOpts.dm_mult)
                    setappdata(0, 'create_mask_gica', 1);
                    return;
                end
            end
        end
    else
        if (strcmpi(modalityType, 'fmri'))
            if (oldOpts.remove_constant_voxels ~= newOpts.remove_constant_voxels)
                setappdata(0, 'create_mask_gica', 1);
                return;
            end
        end
    end
end