function [fileIndex] = icatb_returnFileIndex(numberFiles)
% function that returns the file index naming

% check index
if numberFiles < 10
    fileIndex = ['00', num2str(numberFiles)];
elseif numberFiles < 100
    fileIndex = ['0', num2str(numberFiles)];
else
    fileIndex = num2str(numberFiles);
end
% end for checking
