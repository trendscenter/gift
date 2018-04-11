function files = icatb_fullFile(varargin)
%% Form full file path using the input directory information and the size of
% the files

inputDir = pwd;
files = [];

if mod(nargin, 2) ~= 0
    error('arguments must be passed in pairs');
end

for ii = 1:2:nargin
    if strcmpi(varargin{ii}, 'files')
        files = varargin{ii + 1};
    elseif strcmpi(varargin{ii}, 'directory')
        inputDir = varargin{ii + 1};
    end
end

if ~strcmp(inputDir(end), filesep)
    inputDir = [inputDir, filesep];
end

if ~isempty(files)
    files = [repmat(inputDir, size(files, 1), 1), files];
else
    files = [];
end

