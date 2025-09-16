function [numElectrodes] = icatb_get_num_electrodes(files)
% Count for number of electrodes for each data-set

if isstruct(files)
    temp = files;
    clear files;
    files = str2mat(temp.name);
    clear temp;
end

% Initialise number of electrodes for each data-set
numElectrodes = zeros(1, size(files, 1));

% Loop over files
for nF = 1:size(files, 1)
    currentFile = deblank(files(nF, :));
    data = icatb_loadData(currentFile);
    dims = [size(data, 1), size(data, 2), size(data, 3)];

    if (nF == 1)
        oldDims = dims;
    end

    if (length(find((oldDims == dims) > 0)) ~= length(dims))
        error('Error:DataDim', ['Data dimensions ([%s]) must be the same for all data-sets.', ...
            '\nPlease check the data dimensions for file %s'], ...
            num2str(oldDims), currentFile);
    end
    numElectrodes(nF) = size(data, 4);
end
% End loop over files