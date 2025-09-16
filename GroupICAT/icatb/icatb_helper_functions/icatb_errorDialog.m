function icatb_errorDialog(errorStr, titleStr, windowStyle)
% error dialog box with the specified error string
% and title

if ~exist('titleStr', 'var')
    titleStr = 'Error Found';
end

% window style
if ~exist('windowStyle', 'var')
    windowStyle = 'modal';
end

H = errordlg(errorStr, titleStr); % open error dialog box
set(H, 'windowstyle', windowStyle); % set window style modal
waitfor(H);
