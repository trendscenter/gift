function [varOut] = icatb_check_char(currentVar)
% check whether the variable is a valid character or not
% If it is not a valid character for naming files then return empty

FindChar = find(currentVar == '/' | currentVar == '\' | currentVar == ':' ...
    | currentVar == '*' | currentVar == '?' | currentVar == '"' ...
    | currentVar == '<' | currentVar == '>');

if isempty(FindChar)
    varOut = currentVar;
else
    varOut = [];
end