function icatb_openHelp
%% HTML Help for GIFT
%

% get full path for toolbox
toolboxPath = which('icatb_openHelp.m');

% Get the directory where toolbox is running
[pathstr] = fileparts(toolboxPath);

% Source folder where files are located
sourcePath = fullfile(pathstr, 'icatb_helpManual.htm');

try
    helpHandle = helpdlg('Opening Group ICA Manual...');

    warning off all;
    s = web(sourcePath, '-browser');
    warning on;
    if s ~= 0
        status = system(['firefox ', sourcePath, ' &']);
        if status == 1
            web(sourcePath, '-new');
        else
            close(helpHandle);
            return;
        end
    end

    close(helpHandle);
catch
    disp(lasterr);
end