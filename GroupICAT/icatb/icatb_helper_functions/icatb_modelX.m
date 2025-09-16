function [X] = icatb_modelX(model, numTimePoints, detrendNumber)
% Form regressor matrix
% 
% Inputs:
% 1. model - Model
% 2. numTimePoints - Number of time points
% 3. detrendNumber - Detrend Number
% 
% Output:
% X - Regressor matrix


% Sine and cosine
unitPoint_sin2phi = (2*pi/(numTimePoints - 1));
unitPoint_sin4phi = (4*pi/(numTimePoints - 1));

vec2phi = (0:unitPoint_sin2phi:2*pi)';
vec4phi = (0:unitPoint_sin4phi:4*pi)';

% Ramp function
rampFun = icatb_unitRamp((1:numTimePoints)');
rampFun = detrend(rampFun, 0); % remove the mean before using ramp Function

% Check detrend number
if detrendNumber == 0
    X = [model, ones(numTimePoints, 1)];
elseif detrendNumber == 1
    X = [model, rampFun, ones(numTimePoints, 1)];
elseif detrendNumber == 2
    % steps in sin and cosine
    X = [model, detrend(sin(vec2phi), 0), detrend(cos(vec2phi), 0), rampFun, ...
        ones(numTimePoints, 1)];
else
    X = [model, detrend(sin(vec4phi),0), detrend(cos(vec4phi),0), ...
        detrend(sin(vec2phi),0), detrend(cos(vec2phi),0), ...
        rampFun, ones(numTimePoints, 1)];
end
% End for checking detrend number

