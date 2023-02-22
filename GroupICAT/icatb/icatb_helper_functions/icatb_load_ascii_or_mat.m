function [data] = icatb_load_ascii_or_mat(P)
% Load data of ascii or MAT extension


if size(P, 1) == 1
    P = deblank(P(1, :));
    temp = load(P);

    if isstruct(temp)
        % MAT file
        fieldNames = fieldnames(temp);
        data = getfield(temp, fieldNames{1});
    else
        % Ascii file
        data = temp;
    end

else

    % Loop over files
    for nn = 1:size(P, 1)
        currentFile = deblank(P(nn, :));
        temp = load(currentFile);

        if isstruct(temp)
            % MAT file
            fieldNames = fieldnames(temp);
            tempData = getfield(temp, fieldNames{1});
        else
            % Ascii file
            tempData = temp;
        end
        clear temp;

        if nn == 1
            data = zeros([size(tempData), size(P, 1)]);
        end

        if length(size(tempData)) == 1
            data(:, nn) = tempData;
        else
            data(:, :, nn) = tempData;
        end

    end
    % End loop over files

    data = squeeze(data);

end