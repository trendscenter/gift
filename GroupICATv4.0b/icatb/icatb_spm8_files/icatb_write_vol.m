function icatb_write_vol(V, data)
% Writes the images using the appropriate SPM functions

if isfield(V(1), 'dt')
    V(1).dt(1) = 16;
end

% Add appropriate path to the file
fName = V(1).fname;

V(1).fname = fName;

p = fileparts(V(1).fname);
if (~isempty(p))
    if (exist(p, 'dir') ~= 7)
        [parDir, childDir] = fileparts(p);
        mkdir(parDir, childDir);
    end
end

icatb_spm_write_vol(V, data);

