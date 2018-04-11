function icatb_exit
% Exit group ica

try
    
    if isappdata(0, 'group_ica_modality')
        rmappdata(0, 'group_ica_modality');
    end
    
    % close matlab server
    icatb_closeMatlabServer;
    
    if (exist('gift_finish.m', 'file') == 2) | (exist('gift_finish.m', 'file') == 6)
        fprintf( 'Executing gift_finish ...\n' );
        % execute gift finish file
        gift_finish;
    else
        delete(get(0, 'children'));
    end
    % end for checking finish file
catch
    icatb_displayErrorMsg;
end
