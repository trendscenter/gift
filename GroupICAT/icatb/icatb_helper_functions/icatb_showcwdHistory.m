function [previousDir] = icatb_showcwdHistory

% get the recent previous directories. previousDir contains pwd, GIFT path, and the
% other recent directories

% Pull the previous working directories in Matlab by reading the file
% cwdhistory.m

icatb_defaults;
global STORE_DIRECTORY_INFORMATION;
global DIRS_TO_BE_STORED;

try
    
    % get the old directory
    oldDir = deblank(pwd);
    
    %%%%%%%%%%%%%%%%%%% GIFT Path %%%%%%%%%%%%%%%%%%%%%%%%
    % get the path for the GIFT Toolbox
    whichGift = which('gift');
    giftPath = fileparts(whichGift);
    giftPath = deblank(giftPath);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Get the preferences directory
    preferenecesDir = prefdir;
    
    preferenecesDir = deblank(preferenecesDir);
    
    % Contents of the folder in which preferences are present
    [fileContents] = icatb_listFiles_inDir(preferenecesDir, '*.m');
    
    checkCwdHistory = strmatch('cwdhistory.m', lower(fileContents), 'exact');
    
    % If the cwdhistory file is not present
    % return current working directory and the GIFT Path
    if isempty(checkCwdHistory)
        previousDir{1} = oldDir;
        previousDir{2} = giftPath;
    else
        % open current working directory history m file and read the contents
        fid = fopen(fullfile(preferenecesDir, 'cwdhistory.m'), 'r');
        temp = {};
        if fid ~= -1
            % Read the contents of the cwdhistory and
            [temp, nLines, nonSpace] = icatb_readContents_file(fid);
            fclose(fid);
        end
        
        if isempty(temp)
            error('No directories exist');
        end
        
        
        % Pull only five recent working directories
        if nonSpace > 5
            for ii = 1:5
                previousDir{ii} =  temp{ii};
            end
            % Check if the GIFT path ia already present
            checkGiftPath = strmatch(lower(giftPath), str2mat(lower(previousDir)), 'exact');
            % If GIFT path is not present
            if isempty(checkGiftPath)
                previousDir{6} = giftPath;
            end
        else
            previousDir = temp;
            if isempty(checkGiftPath)
                previousDir{length(temp) + 1} = giftPath;
            end
        end
    end
    
    % preserve the old directory
    cd(oldDir);
    
    % Check the present working directory in the previous directories
    checkPwd = strmatch(lower(oldDir), str2mat(lower(previousDir)), 'exact');
    
    % Add present working directory
    if isempty(checkPwd)
        previousDir{length(previousDir) + 1} = oldDir;
    end
    
    count = length(previousDir);
    % Subject directories can be stored
    if strcmp(lower(STORE_DIRECTORY_INFORMATION), 'yes')
        for ii = 1:length(DIRS_TO_BE_STORED)
            if ~isempty(DIRS_TO_BE_STORED{ii})
                count = count + 1;
                previousDir{count} = DIRS_TO_BE_STORED{ii};
            end
        end
    end
    
catch
    
    previousDir{1} = oldDir;
    previousDir{2} = giftPath;
    count = length(previousDir);
    % Subject directories can be stored
    if strcmp(lower(STORE_DIRECTORY_INFORMATION), 'yes')
        for ii = 1:length(DIRS_TO_BE_STORED)
            if ~isempty(DIRS_TO_BE_STORED{ii})
                count = count + 1;
                previousDir{count} = DIRS_TO_BE_STORED{ii};
            end
        end
    end
    
    % preserve the old directory
    cd(oldDir);
end
