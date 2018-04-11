function varOut = icatb_constructString(varIn)
% construct number to string 

minVar = min(varIn);
maxVar = max(varIn);
% max of variable is 0 and min of variable is 0
if (minVar == maxVar)
    varOut = num2str(minVar);
else
    varOut = [num2str(varIn(1)), ':', num2str(mean(diff(varIn))), ':', num2str(varIn(end))];
end