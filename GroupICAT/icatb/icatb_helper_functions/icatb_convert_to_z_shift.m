function icatb_convert_to_z_shift(files)
% Function to convert images to z-shift

% Load defaults
icatb_defaults;
global SMOOTHINGVALUE;

if ~exist('files', 'var')
    files = icatb_selectEntry('typeSelection', 'multiple', 'typeEntity', 'file', 'fileType', 'image', 'title', 'Select images ...', 'filter', ...
        '*.img;*.nii', 'filenumbers', 1);
end

if isempty(files)
    error('Files are not selected to center image/images distribution to zero');
end

drawnow;

files = str2mat(files);

% Rename files by adding a file number at the end
files = icatb_rename_4d_file(files);

disp('Centering image/images distribution to zero ...');

% Loop over files
for nF = 1:size(files, 1)

    % Current file
    file_name = deblank(files(nF, :));

    [pathstr, name, extn] = fileparts(file_name);

    if isempty(pathstr)
        pathstr = pwd;
    end

    [extn, number] = icatb_parseExtn(extn);

    V = icatb_spm_vol(file_name);

    %     fileIndex = '';
    %
    %     if (number > 1)
    %         fileIndex = ['_', icatb_returnFileIndex(number)];
    %     end

    % Read data
    data = icatb_spm_read_vols(V);

    % Output file info
    Vout = V(1);
    Vout.fname = fullfile(pathstr, ['Zs_', name, extn]);
    Vout.n(1) = number;

    % Convert data to z-score
    data = data / stdN(data(data ~= 0));

    % Bins
    hist_ind = min(data(:)):0.1:max(data(:));

    % Compute histogram
    data_hist = hist(data(data ~= 0), hist_ind);

    data_hist = data_hist(:);

    [Y, I] = max(icatb_gauss_smooth1D(data_hist, SMOOTHINGVALUE));

    data = (data ~= 0).*(data - hist_ind(I));

    data = data / std(data(data ~= 0));

    % Write Z-shift
    icatb_write_vol(Vout, data);

end
% End of loop over files

disp('Done');

fprintf('\n');

disp('New images are named with prefix Zs_');

fprintf('\n');


function out = stdN(in)
%std of image or matrix

out = std(in(:));