function icatb_checkVariable(inputData, keywd, file_name)
% error check for the required field in the structure

if isempty(getfield(inputData, keywd))
    error(['Undefined variable or not the required data type (', keywd, ') in file: ', file_name]);
end
