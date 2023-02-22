function [newVar] = icatb_check_variable_integer(currentVar)
 
 % check whether the variable in the string is integer or not

 % check characters not caught by the str2num and the floating point
 checkNumber = find(currentVar == '''' | currentVar == ',' | currentVar == ';' | ...
                    currentVar == 'i' | currentVar == 'j');
 
                
% check if the string is logical
if strcmp(lower(currentVar), 'true') | strcmp(lower(currentVar), 'false')
    checkLogical = currentVar;
else
    checkLogical = [];
end

% Vars should not contai logical, floating point format
% str2num takes care of special characters and alphabets
 if isempty(checkNumber) & isempty(checkLogical)
     newVar = str2num(currentVar);     
     diff_var = newVar - round(newVar);
     if(diff_var ~= 0)
         newVar = [];
     end
 else
     newVar = [];
 end