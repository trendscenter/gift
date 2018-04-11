function noisecloud_spm_coregister(coregRefImage, coregSourceImage, otherImages, outDir)
%% Coregister and reslice given images to the reference image
%
% Inputs:
% coregRefImage - Reference image
% coregSourceImage - Source image
% otherImages - Any other images that need to be transformed to reference
% image space
% outDir - Output directory of resliced images
%
% Output images are written to the disk with r* prefix.
%

coreg = 1;
reslice = 1;

verFlag = lower(spm('ver'));
verNum = str2num(strrep(verFlag, 'spm', ''));

if (verNum >= 8)
    % spm 8 and 12
    defaults = spm_get_defaults;
else
    % spm 5
    spm_defaults;
    global defaults;
end

if (~exist('otherImages', 'var'))
    otherImages = '';
end

eoptions = defaults.coreg.estimate;

%% Coregister
% The above code is from spm_config_coreg function
disp(['Coregistering ', coregSourceImage, ' to ', coregRefImage, ' ...']);
x  = spm_coreg(coregRefImage, coregSourceImage, eoptions);
M  = inv(spm_matrix(x));

%% Get orientation
PO = strvcat(coregSourceImage, otherImages);
MM = zeros(4,4,size(PO,1));
for j = 1:size(PO,1),
    MM(:, :, j) = spm_get_space(deblank(PO(j, :)));
end

%% Apply transformation to given images
files_to_reslice = cell(1, size(PO, 1));
header_files =  cell(1, size(PO, 1));
for j = 1:size(PO, 1)
    fname = deblank(PO(j, :));
    [pDir, fN, extn] = fileparts(fname);
    if (~exist('outDir', 'var'))
        outDir = pDir;
    end
    V = spm_vol(fname);
    data = spm_read_vols(V);
    V.mat = M*MM(:, :, j);
    commaPos = find(extn == ',');
    if (~isempty(commaPos))
        extn = extn(1:(commaPos-1));
    end
    V.fname = fullfile(outDir, [fN, '_', extn]);
    files_to_reslice{j} = V.fname;
    spm_write_vol(V, data);
    if (strcmpi(extn, '.img'))
        header_files{j} = fullfile(outDir, [fN, '_', '.hdr']);
    end
    % spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
end;

files_to_reslice = unique(files_to_reslice);
if (strcmpi(extn, '.img'))
    header_files = unique(header_files);
end

chk = cellfun('isempty', header_files);
header_files(chk) = [];

%% Reslice
disp('Reslicing ...');
P = strvcat(coregRefImage, char(files_to_reslice));
% Reslice options used from spm_defaults
flags.mask   = defaults.coreg.write.mask;
flags.mean   = 0;
%flags.interp = defaults.coreg.write.interp;
flags.interp=0;
flags.which  = 1;
flags.wrap   = defaults.coreg.write.wrap;

spm_reslice(P, flags);

%% Cleanup intermediate images and headers
for nF = 1:length(files_to_reslice)
    delete (files_to_reslice{nF});
end

if (strcmpi(extn, '.img'))
    for nF = 1:length(header_files)
        delete (header_files{nF});
    end
end

%% Remove underscore in resliced file names
files = strvcat(char(files_to_reslice), char(header_files));
for nF = 1:size(files, 1)
    fname = deblank(files(nF, :));
    [p, fN, extn] = fileparts(fname);
    oF = fullfile(p, ['r', fN, extn]);
    fN(end) = []; % remove underscore
    cF = fullfile(p, ['r', fN, extn]);
    movefile(oF, cF);
end

fprintf('Done\n');