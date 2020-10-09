function inputVector = icatb_unitRamp(inputVector)

% calculates the unit ramp function
maxVal = max(inputVector);

inputVector = inputVector/maxVal;