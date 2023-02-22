function icatb_openHTMLHelpFile(fileName)
%% Opens the specified HTML help file in the browser
%

setupICAFilePath = which(fileName); % get the full file name
if (isempty(setupICAFilePath))
    setupICAFilePath = fileName;
end
[pathstr, fileName, extn] = fileparts(setupICAFilePath); % get the path for the html file
fullFileName = fullfile(pathstr, [fileName, extn]); % form the full file path

warning off all;
s = web(fullFileName, '-browser'); % open HTML file in browser
warning on;

if s ~= 0
    status = system(['firefox ', fullFileName, ' &']);
    if status == 1
        s = web(fullFileName, '-new');
    end
end