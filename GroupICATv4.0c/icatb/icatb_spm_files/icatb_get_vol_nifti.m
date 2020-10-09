function V = icatb_get_vol_nifti(file)
% get volume of file by parsing the extension

% Parse extension
[file, fileNum] = icatb_parseExtn(file);

V = icatb_spm_vol_nifti(file, fileNum);