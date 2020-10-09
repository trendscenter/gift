function icatb_closeMatlabServer
% close Matlab server on PC's for nojvm

global COMSERVERHANDLE;

try
    % Quit COM server only for windows
    if ispc
        if ~usejava('jvm')
            disp('Closing Matlab Server ...');
            fprintf('\n');
            Quit(COMSERVERHANDLE);
            delete(COMSERVERHANDLE);
            clear global COMSERVERHANDLE;
        end
    end

catch

end