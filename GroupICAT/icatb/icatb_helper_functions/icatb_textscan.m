function data = icatb_textscan(file_name)
%% Read text file.
%
% Inputs:
% 1. file_name - File name
%
% Outputs:
% data - cell array of strings
%

% Open file in read mode
fid = fopen(file_name, 'r');
if (fid == -1)
    error('Error:txtfile', 'File %s cannot be opened\n', file_name);
end

tempByte = fread(fid, [1, 1], 'char');
frewind(fid);

if (tempByte >= 255)
    fprintf('\n');
    disp(['File ', file_name, ' is in unicode format. Bytemarkers will be excluded while reading.']);
    fprintf('\n');
    fclose(fid);
    warning off MATLAB:iofun:UnsupportedEncoding;
    fid = fopen(file_name, 'r', 'l', 'UTF16-LE');
    fseek(fid, 2, 0);
    % Read all strings in one cell
    data = textscan(fread(fid, '*char'), '%s', 'delimiter', '\n', 'whitespace', '');
else
    % Read all strings in one cell
    data = textscan(fid, '%s', 'delimiter', '\n', 'whitespace', '');
end

fclose(fid);

data = data{1};