function [timecourse, topography] = icatb_eeg_loadICAData(varargin)
% load EEG components

% TC reshape
tcReshape = 'yes';
moving_average_window = 3;

% get the required parameters
for ii = 1:2:nargin
    % Check parameters
    if strcmpi(varargin{ii}, 'structdata')
        % structural file
        structData = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'compfiles')
        % component files
        compFiles = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'comp_numbers')
        % component numbers
        comp_numbers = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'tcreshape')
        % reshape tc to 3D
        tcReshape = varargin{ii + 1};
    end
    % end for checking parameters
end
% end for loop

if ~exist('compFiles', 'var')
    error('Component files are missing');
end

% Structural dimensions
structDIM = [size(structData, 1), size(structData, 2), size(structData, 3)];

load(compFiles);

timecourse = timecourse(comp_numbers, :);
topography = topography(:, comp_numbers);

chkTrials = (structDIM(2) >= (2*moving_average_window + 1));

if (chkTrials)
    if strcmpi(tcReshape, 'yes')
        timecourse = reshape(timecourse, [size(timecourse, 1), structDIM]);
        
        % Loop over components
        for nComp = 1:size(timecourse, 1)
            % Loop over samples
            for nSamp = 1:size(timecourse, 2)
                timecourse(nComp, nSamp, :) = icatb_moving_average(timecourse(nComp, nSamp, :), moving_average_window);
            end
            % End loop over samples
        end
        % End loop over components
        
    end
end