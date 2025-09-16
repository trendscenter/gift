function icatb_displayErrorMsg
% Display error message and print the stack information on Matlab 7 and
% higher

tempErr = lasterror; % Get the last error

% If stack field is present print stack info
if isfield(tempErr, 'stack')
    fprintf('\n');  
    disp('......................................');
    msgString(1).string = sprintf('%s\n', 'Group ICA Error Information: ');
    msgString(length(msgString) + 1).string = sprintf('%s\n', tempErr.message); % Error message
    % Print the line numbers
    for ii = 1:length(tempErr.stack)
        msgString(length(msgString) + 1).string = sprintf('Error in ==> <a href="matlab:opentoline(''%s'', %d)">%s at %d</a>', ...
            tempErr.stack(ii).file, tempErr.stack(ii).line, tempErr.stack(ii).name, ...
            tempErr.stack(ii).line);
    end
    % display the message and produce beep sound
    disp(str2mat(msgString.string));
    disp('......................................');
    clear tempErr;
    error(' ');
else
    rethrow(lasterror);
end