function icatb_openMatlabServer
% open Matlab server on PC's for nojvm

global COMSERVERHANDLE;

if ispc
    % no java
    if ~usejava('jvm')
        % check Matlab server handle
        if ~(checkCOM(COMSERVERHANDLE))
            % open matlab COM automation server
            comOpenStr = 'Opening Matlab server ...';
            disp(comOpenStr);
            comHelpH = helpdlg(comOpenStr, 'Matlab Server');
            COMSERVERHANDLE = actxserver('matlab.application');
            set(COMSERVERHANDLE, 'visible', 0); % make it not visible
            try
                delete(comHelpH);
            catch
            end
        end
        % check if it is Matlab server handle
        
    end
    % no java
end
% for pc

function status = checkCOM(h)
% Detect if it is COM object

status = 0;

if ~isempty(h)
    versionNum = icatb_get_matlab_version;
    
    if versionNum < 14
        status = isa(h, 'COM.matlab.application');
    else
        status = iscom(h);
    end
end