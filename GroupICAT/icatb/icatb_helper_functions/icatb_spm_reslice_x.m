function s_spatial_ref = icatb_spm_reslice_x(firstFile, s_source_template, s_workdir)
    % Wrapper function that adds little GIFT extra for GIFT dialekt of spm_reslice
    % Cyrus Eierud
    % Copy template to workdir
    s_source_template = icatb_spm_file(s_source_template, 'number', ''); %strip vol num
    [pathstr, fp, extn] = fileparts(s_source_template);
    fileContents = icatb_listFiles_inDir(pathstr, [fp, extn]);
    s_file_nmark_del = [s_workdir filesep fileContents];
    copyfile(s_source_template, s_file_nmark_del)

    % Create list of volumes to send into icatb_spm_reslice
    V_tmp = spm_vol(s_file_nmark_del);
    nVols = numel(V_tmp);
    cs_nmark_list = cell(nVols+1,1);
    cs_nmark_list{1} = firstFile;
    for i = 2:nVols+1
        cs_nmark_list{i} = sprintf('%s,%d', s_file_nmark_del, i-1);
    end

    % Execute icatb_spm_reslice
    flags.mask   = 0;
    flags.mean   = 0;
    flags.interp = 4;
    flags.which  = 1;
    flags.prefix  = 'r';
    flags.wrap   = [1, 1, 0];
    icatb_spm_reslice(cs_nmark_list,flags);
    delete(s_file_nmark_del);
    s_spatial_ref = [s_workdir filesep flags.prefix fileContents];
end

