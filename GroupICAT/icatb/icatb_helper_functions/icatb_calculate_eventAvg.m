function eventAvg = icatb_calculate_eventAvg(tempTC, interpFactor, TR, windowSize, selectedOnset)
% calculates event average based on the onset, window size, interpolation
% factor
%
% Inputs:
% 1. tempTC - Time course
% 2. interpFactor - Interpolation factor
% 3. TR 
% 4. windowSize (sec)
% 5. selectedOnset - vector of onsets
%
% Output:
% eventAvg - event average

% Original Interpolation factor
origInterpFactor = interpFactor;

% Interpolation factor
interpFactor = ceil(origInterpFactor*TR);

% Interpolated ICA
interpICA = icatb_interp(tempTC, interpFactor);

eventAvg = zeros(windowSize*origInterpFactor, 1); % initialise event average

%% Event Related Average
for ons = 1:length(selectedOnset)
    if selectedOnset(ons)*interpFactor == 0
        startOnset = 1;
    else
        startOnset = selectedOnset(ons)*interpFactor; % start onset
    end
    endOnset = startOnset + windowSize*origInterpFactor - 1; % end onset
    if endOnset <= length(interpICA)
        eventAvg = eventAvg + interpICA(startOnset:endOnset);
        lastNumber = ons;
    end
end

if exist('lastNumber', 'var')
    eventAvg = eventAvg / lastNumber; % event Average
else
    eventAvg = ones(size(eventAvg))*interpICA(startOnset);
end