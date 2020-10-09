function [icasig, A] = icatb_eeg_scaleComp_dataUnits(origData, icasig, A, DIM)
% Convert eeg components to data units by fitting average component signal
% to average EEG signal.
%
% Input:
% 1. origData - Original data is a 2D matrix of dimensions (time*trials by electrodes) by
% electrodes
% 2. icasig - Component data is a 2D matrix of dimensions (time by trials) by
% components
% 3. A - Mixing matrix is a 2D matrix of dimensions components by
% electrodes
% 4. DIM - Dimensions of data
%
% Output:
% 1. icasig - Scaled component data
% 2. A - Scaled mixing matrix

icatb_defaults;

global DETRENDNUMBER;

A = A';


if (prod(DIM(1:2)) == size(origData, 1))
    % Reshape data as time by trials by electrodes
    origData = reshape(origData, [DIM(1:2), size(origData, 2)]);
end


%%%% Scale components using original data as observation %%%%
for nComp = 1:size(A, 2)
    % Current component
    compData = icasig(nComp, :);
    if prod(DIM(1:2)) == length(compData)
        compData = reshape(compData, DIM(1:2));
        % Average component signal
        compData = mean(compData, 2);
        % Find the electrode that contributes maximum weight to the current
        % component
        [dd, ind] = max(abs(A(:, nComp)));
        tempData = squeeze(origData(:, :, ind));
        % Average signal at the maximum electrode
        tempData = mean(tempData, 2);
    else
        % Find the electrode that contributes maximum weight to the current
        % component
        [dd, ind] = max(abs(A(:, nComp)));
        tempData = origData(:, ind);
    end

    % Compute regression using data as observation and component as model
    [rSquare_stat, b, ModelIndices, otherIndices] = icatb_multipleRegression(compData, tempData, 1, 1, length(tempData), DETRENDNUMBER);
    b = b(1);
    b(isnan(b)) = 0;

    if abs(b) < 1e-5
        warning(['Regression slope (beta weight) of component ', num2str(nComp), ' while converting to data units is close to zero']);
    end

    % Scale components based on regression slope
    icasig(nComp, :) = abs(b)*icasig(nComp, :);
    A(:, nComp) = abs(b)*A(:, nComp);

end
%%%% End for scaling components using original data as observation %%%%