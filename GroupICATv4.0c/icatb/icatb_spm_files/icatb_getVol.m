function V = icatb_getVol(files, number)
%% Get volume
%

files = icatb_rename_4d_file(files);

if (exist('number', 'var'))
    if (max(number) > size(files, 1))
        error('Image number requested exceeds the number of images');
    end
    files = deblank(files(number, :));
end

V = icatb_spm_vol(files);
