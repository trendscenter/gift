function [pos] = icatb_findstr(str, pattern)
% finds the string based on the pattern

check1 = which('findstr.m');

check2 = which('strfind.m');

if isempty(check1)
    % if there is no find string
    pos = strfind(str, pattern);
else
    pos = findstr(str, pattern);
end