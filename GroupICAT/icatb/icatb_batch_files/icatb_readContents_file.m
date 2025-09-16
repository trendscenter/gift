function [readContents, count, nonSpace] = icatb_readContents_file(fid)

% Reads the contents of the files
% Reads empty characters

% Output: 
%        1. readContents: text of the file not containing empty spaces in
%        new line
%        2. Number of lines of  the file

count = 0;

% loop until end of file is encountered
while ~feof(fid)       
    % Read the current line
    tline = fgets(fid);  
    tline = deblank(tline);
    count = count + 1;
    temp{count} = tline;                       
end

% Place the position of the file at the beginning of the record
frewind(fid);

nonSpace = 0;
for ii = 1:count
    if ~strcmp(temp{ii}, '')
        nonSpace = nonSpace + 1;
        readContents{nonSpace} = deblank(temp{ii});
    end
end