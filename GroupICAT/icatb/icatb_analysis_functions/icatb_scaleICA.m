function [varargout] = icatb_scaleICA(icasig, A, fmriFiles, scaleType, dataType, complexInfo, mask_ind)
% [scaledICASig scaledA]=scaleICA(icasig,A,origData)
% ---------------------------------------------------------------
% This function takes the timecourses and the images resulting from ICA. It also
% takes the original data that was Analyzed with ICA.
% The function returns a scaled timecourse and a scaled image.
% -----------------------------------------------------------------------
% -icasig = ica signal in 2D(number of components, x*y*z)
% -A = ica timecourse in 2d(number of components, number of timepoints)
% -origData = original dataata that icasig and A was obtained from.
% data can be either 2D(x*y*z,number of timepoints) or 4D(x,y,z,t)
% -scaleType lets you select how to scale the timecourse
% 1 means difference between max and min
% 2 means the top 10% and bottom 10%
% fmriFiles = functional files to load
% --------------------------------------------------------------------

% load defaults
icatb_defaults;
global WRITE_COMPLEX_IMAGES;
global SCALE_DEFAULT;
global SBM_Z_SCORES_RMMEAN;

if (isempty(SBM_Z_SCORES_RMMEAN))
    SBM_Z_SCORES_RMMEAN = 1;
end

[modalityType, dataTitle] = icatb_get_modality;

% Calibrate Options
if (strcmpi(modalityType, 'fmri') || strcmpi(modalityType, 'smri'))
    calibrateOptions = char('No Scaling', 'Scale to Original Data(%)', 'Z-scores', 'Scaling in Timecourses', 'Scaling in Maps and Timecourses');
else
    calibrateOptions = char('No Scaling', 'Scale to Original Data', 'Z-scores', 'Scaling in Topographies', 'Scaling in Timecourses and Topographies');
end

% Add one to the index
if ~exist('scaleType', 'var')
    scaleType = SCALE_DEFAULT;
end

scaleType = scaleType + 1;

if scaleType > size(calibrateOptions, 1)
    error(['Presently there are ', num2str(size(calibrateOptions, 1)), ' scaling methods only.']);
end

% calibrate Type
calibrateType = lower(deblank(calibrateOptions(scaleType, :)));

% Return the available calibration schemes if no input arguments are
% specified
if nargin == 0
    varargout{1} = calibrateOptions;
    return;
end

if ~exist('dataType', 'var')
    dataType = 'real';
end

if ~exist('complexInfo', 'var')
    complexInfo = [];
end

if strcmpi(dataType, 'complex')
    % convert component images to complex data
    icasig = complex_data(icasig);
    % convert time courses to complex data
    A = complex_data(A);
end

if (strcmpi(calibrateType, 'no') || strcmpi(calibrateType, 'no scaling'))
    % Not calibrating
    A = A'; % Transposing A such that it is of dimension Timepoints by components (each column represents a time course)
    disp('Not scaling components');

elseif (strcmpi(calibrateType, 'calibrate') || strcmpi(calibrateType, 'scale to original data(%)') || strcmpi(calibrateType, 'scale to original data'))
    % calibrating components to percent signal change

    %disp('--Loading Functional Data and Calibrating Components');
    disp(['--Loading ', dataTitle, ' data and scaling components to data units']);

    % load original data for scaling timecourses and
    % components
    [origData, HInfo] = icatb_loadData(fmriFiles, dataType, complexInfo, 'read', [], mask_ind);

    %[origData, HInfo] = icatb_read_data(fmriFiles, [], mask_ind);

    % Adjust global variable to make components as complex_data type
    if strcmpi(dataType, 'complex')
        origWriteImages = WRITE_COMPLEX_IMAGES;
        WRITE_COMPLEX_IMAGES = complexInfo.complexType;
        % convert original data to complex data
        origData = complex_data(origData);
        % restore the original value for the WRITE_COMPLEX_IMAGES
        WRITE_COMPLEX_IMAGES = origWriteImages;
    end

    %force timecourse to be components*voxels
    if size(A,1) > size(A,2)
        A=A';
    end

    if strcmpi(modalityType, 'fmri')
        [icasig, A] = icatb_percentSignalChange(origData, icasig, A);
    else
        % Scaled component data
        [icasig, A] = icatb_eeg_scaleComp_dataUnits(origData, icasig, A, HInfo.DIM);
    end

elseif (strcmpi(calibrateType, 'z-scores'))
    % calibrating to Z scores both the components and timecourses

    disp('Converting components to z-scores');

    if strcmpi(modalityType, 'fmri')
        A = icatb_convertToZScores(A); % time course
        icasig = icatb_convertImageToZScores(icasig); % spatial maps
    elseif strcmpi(modalityType, 'smri')
        % Loop over components
        for n = 1:size(A, 1)
            if (SBM_Z_SCORES_RMMEAN)
                 A(n, :) = detrend(A(n, :), 0);
            end
            A(n, :) = A(n, :)./std(A(n, :));
        end
        % End loop over components
        icasig = icatb_convertImageToZScores(icasig); % spatial maps
    else

        % Loop over components
        for n = 1:size(A, 1)
            A(n, :) = A(n, :)./std(A(n, :));
            icasig(n, :) = icasig(n, :)./std(icasig(n, :));
        end
        % End loop over components
    end

    A = A';

elseif (strcmpi(calibrateType, 'scaling in timecourses') || strcmpi(calibrateType, 'scaling in topographies') || ~isempty(findstr(calibrateType, 'norm ics')))
    % Put max intensity information in independent components

    if (~strcmpi(modalityType, 'eeg'))
        disp('Spatial maps are normalized using the average of top 1% voxels and the resultant value is multiplied to the timecourses ...');
    else
        disp('Timecourses are normalized using the average of top 1% amplitudes and the resultant value is multiplied to the topographies ...');
    end

    for n = 1:size(icasig, 1)
        %maxVal = max(icasig(n, :));
        tmp = icasig(n, :);
        tmp = tmp(tmp > 0);
        maxVal = 1;
        if (~isempty(tmp))
            voxels = ceil(0.01*length(icasig(n, :)));
            voxels = min([voxels, length(tmp)]);
            [dd, inds] = sort(tmp);
            maxVal = mean(tmp(inds(end:-1:end-voxels+1)));
        end
        A(n, :) = A(n, :)*maxVal;
        icasig(n, :) = icasig(n, :)/maxVal;

    end
    A = A';

elseif (strcmpi(calibrateType, 'scaling in maps and timecourses') || strcmpi(calibrateType, 'scaling in timecourses and topographies') || ...
        ~isempty(findstr(calibrateType, 'joint scaling')))
    % Put combined scaling information in both independent components and
    % timecourses

    if (~strcmpi(modalityType, 'eeg'))
        disp('Spatial maps are scaled using the standard deviation of timecourses and timecourses are scaled using the maximum spatial intensity value ...');
    else
        disp('Timecourses are scaled using the standard deviation of topographies and topographies are scaled using the maximum amplitudes of timecourses ...');
    end

    for n = 1:size(icasig, 1)
        maxVal = max(icasig(n, :));
        stdVal = std(A(n, :));
        A(n, :) = A(n, :)*maxVal;
        icasig(n, :) = icasig(n, :)*stdVal;
    end
    A = A';

    % Add your calibrate code below

end
% end for scaling steps

scaledA = A; % scaled time course of dimension Time points by components
scaledICASig = icasig; % scaled component of dimension components by voxels

% return output
varargout{1} = scaledICASig;
varargout{2} = scaledA;