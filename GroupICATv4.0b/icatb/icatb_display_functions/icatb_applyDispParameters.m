function [icasig] = icatb_applyDispParameters(icasig, convertToZ, returnValue, threshValue, ...
    structDIM, HInfo)
% apply display parameters to the component images
% check whether the image is of type structure, complex or real

% check if the image is of class name complex_data
if isa(icasig, 'complex_data')
    
    % get the field names
    fieldNames = fieldnames(icasig);
    
    % set the fields to the object
    for ii = 1:length(fieldNames)
        % get the field data
        currentData = getfield(icasig, fieldNames{ii});
        % apply display parameters
        currentData = icatb_applyDispParameters_comp(currentData, convertToZ, returnValue, threshValue, ...
        structDIM, HInfo);
        % set the data to the object
        icasig = setfield(icasig, fieldNames{ii}, currentData);
    end

elseif isreal(icasig)
    
    % apply display parameters to images
    [icasig] = icatb_applyDispParameters_comp(icasig, convertToZ, returnValue, threshValue, ...
        structDIM, HInfo);

else

    error('Unknown data type');

end
% end for checking the image data type