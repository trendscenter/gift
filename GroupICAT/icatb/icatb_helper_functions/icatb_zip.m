function icatb_zip(zipFile, files_in_zip, targetDir)
% zips the file

global COMSERVERHANDLE;

oldDir = pwd;

if ~exist('targetDir', 'var')
    targetDir = pwd;
end

if ~exist('zipFile', 'var')
    error('No zip file name passed');
end

if ~exist('files_in_zip', 'var')
    error('files_in_zip variable should be passed ');
end


% for java use Matlab function
if usejava('jvm')
    % change to target directory
    cd(targetDir);
    % matlab function for zip
    zip(zipFile, files_in_zip);
    % change to old directory
    cd(oldDir);
else

    % check the OS
    if ispc
        % make adjustment for the nojvm case

        % mat file to store zip file information
        zipMatFile = fullfile(targetDir, 'zipInfo.mat');
        % save the zip file information like zip file name, files in zip and
        % target directory
        save(zipMatFile, 'zipFile', 'files_in_zip', 'targetDir', 'oldDir');
        % put character array in Matlab server
        PutCharArray(COMSERVERHANDLE, 'zipMatFile', 'base', zipMatFile);
        Execute(COMSERVERHANDLE, 'load(zipMatFile); cd(targetDir); zip(zipFile, files_in_zip); delete(zipMatFile); cd(oldDir)');
        % change to old directory
        cd(oldDir);
        if ~exist(zipFile, 'file')
            error(['Failed to create the zip file: ', zipFile]);
        end


    else
        % run system command for other operating systems

        % delete the zip file first
        if exist(zipFile, 'file')
            delete(fullfile(targetDir, zipFile));
        end

        cd(targetDir);

        % command to enter
        commandToEnter = ['zip "', zipFile, '"'];

        % loop over files
        for ii = 1:length(files_in_zip)
            % command to enter
            commandToEnter = [commandToEnter, ' "', files_in_zip{ii}, '"'];
        end
        % end for loop

        % run system command
        [checkStatus, result] = system(commandToEnter);

        % change to old directory
        cd(oldDir);

        % if the command failed
        if checkStatus == 1
            error(result);
        end
        % do error checking

    end
    % end for checking the OS

end
% end for checking java
