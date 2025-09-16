function icatb_unzip(zipFile, targetDir)
% unzip the file

global COMSERVERHANDLE;

oldDir = pwd;

if ~exist('targetDir', 'var')
    targetDir = pwd;
end

if ~exist('zipFile', 'var')
    error('No zip file name passed');
end


% with java
if usejava('jvm')
    cd(targetDir);
    % Matlab function for unzipping
    unzip(zipFile);
    cd(oldDir);
else

    % for windows
    if ispc

        % make adjustment for the nojvm case
        if ~exist(zipFile, 'file')
            error(['zip file: ', zipFile, ' doesn''t exist']);
        end

        % put character array in Matlab server
        PutCharArray(COMSERVERHANDLE,  'zipFile', 'base',  zipFile);
        PutCharArray(COMSERVERHANDLE, 'targetDir', 'base', targetDir);
        PutCharArray(COMSERVERHANDLE, 'oldDir', 'base', oldDir);
        Execute(COMSERVERHANDLE, 'cd(targetDir); unzip(zipFile); cd(oldDir);');
        cd(oldDir);

    else
        % for other operating systems

        cd(targetDir);

        commandToEnter = ['unzip "', zipFile, '"'];
        % run unzip command
        [checkStatus, result] = system(commandToEnter);

        cd(oldDir);

        if checkStatus == 1
            error(result);
        end
        % end for error checking
    end
    % end for checking OS

end
% end for checking java or nojava
