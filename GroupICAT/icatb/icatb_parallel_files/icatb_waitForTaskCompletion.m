function icatb_waitForTaskCompletion(files, currentTime)
%% Waitfor task completion (Synchronize processes)
%

files = cellstr(files);

chkFiles = ones(1, length(files));
filesVec = (1:length(files));
while 1
    for nF = filesVec
        tmp = dir(files{nF});
        if (~isempty(tmp) && (tmp.datenum >= currentTime))
            chkFiles(nF) = 0;
        end
    end
    filesVec = find(chkFiles);
    if (isempty(filesVec))
        break;
    end
end