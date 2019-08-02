function [file, number] = icatb_parseExtn(file)
% Parse extension

[pathstr, fName, extn] = fileparts(file);

commaPos = find(extn == ',');

number = 1;

if ~isempty(commaPos)
    number   = str2num(extn((commaPos(1)+1):end));
    extn = extn(1:(commaPos-1));
    file   = fullfile(pathstr, [deblank(fName), extn]);
end
