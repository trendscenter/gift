function [eventAvg, convolved_tc, selectedOnset] = icatb_calculate_eventAvg_deconv(tempTC, TR, windowSize, selectedOnset)
% Calculates event average using deconvolution funtion written by Tom
% Eichele
%
% Inputs:
% 1. tempTC - Time course
% 2. TR - Response time
% 3. windowSize - Window size
% 4. selectedOnset - vector of onsets
%
% Outputs:
% 1. eventAvg - event average
% 2. convolved_tc - Convolved Timecourse

if (selectedOnset(1) == 0)
    selectedOnset(1) = 1;
end

% Convert to row vector
tempTC = tempTC(:)';

% Hypothetical TC
tc_hyp = zeros(length(tempTC), 1);
tc_hyp(selectedOnset) = 1;
tc_hyp = tc_hyp(:)';

% Convert tc_hyp to z-scores
tc_hyp = detrend(tc_hyp, 0);
tc_hyp = tc_hyp./std(tc_hyp);

% Convert TC to z-scores
tempTC = detrend(tempTC, 0);
tempTC = tempTC./std(tempTC);


tempTC = tempTC(:);
tc_hyp = tc_hyp(:);

% HRF
myHrf = icatb_spm_hrf(TR);

if length(myHrf) > windowSize
    myHrf = myHrf(1:windowSize);
end

% Convolve TC with HRF
convolved_tc = conv(tempTC, myHrf);

% Deconvolve
%eventAvg = icatb_deconvolve(convolved_tc, tc_hyp);

% Deconvolve only the measured time course not the convolved time course
eventAvg = icatb_deconvolve(tempTC, tc_hyp(1:end-windowSize));

eventAvg = eventAvg(1:windowSize);