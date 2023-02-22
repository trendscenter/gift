function StatusHandle = icatb_statusBar(cmd_line,varargin)
% StatusHandle = icatb_statusBar(cmd_line,varargin)
%--------------------------------------------------------------------------
% CREATED BY: Eric Egolf 
% LAST MODIFIED: 12-29-03
% ABOUT: Statusbar figure tellling user how much of a job is done
%
% #########################################################################
% 
% USAGE:
%
% -StatusHandle = icatb_statusBar;
%   INFO: Initialize Generic status bar
%   PARAMETERS: None
%
% -StatusHandle = icatb_statusBar('init',endNumber,title,xlabel,ylabel);
%   INFO: Initialize Status Bar, Only one status bar can be open at once
%   else and error is displayed
%   PARAMETERS:
%   endNumber = number that equals 100% when the status bar counter gets to it
%   title,xlabel,ylabel = descriptive strings
%
% -StatusHandle = icatb_statusBar('increment',incrementAmount);
%  INFO: Increments statusbar
%  PARAMETERS:
%  incrementAmount = amount to increases the counter to
%
% -StatusHandle = icatb_statusBar('finished');
%  INFO: closes status handle
%  PARAMETERS: None
%
% #############################################################
% 
% LICENSING:
% 
%------------------------------------------------------



%get parameters from varargin
if nargin == 0
    cmd_line = 'init';
    endNumber =100;
    title = '';
    xlabel = '';
    ylabel = '';
else
    if(strcmp(cmd_line,'init'))
        endNumber = varargin{1};
        title = varargin{2};
        xlabel = varargin{3};
        ylabel = varargin{4};
    end
    if(strcmp(cmd_line, 'increment'))
        incrementAmount = varargin{1};
    end
    if(strcmp(cmd_line,'finished'));
        %no parameters
    end
end



switch cmd_line
    
    case 'init'
        icatb_defaults;
        global BG_COLOR;
        global FG_COLOR;
        global AXES_COLOR;
        global FONT_COLOR;
        
        
        
        StatusHandle = findobj('Tag','StatusBar');
        if(~isempty(StatusHandle))
            close(StatusHandle);    
        end
        
        statusPos = [ .1 .1 .2 .3];
        visible = 'on';
        StatusHandle = icatb_getGraphics('Status Bar','statusBar','StatusBar');
        set(StatusHandle,'units','normalized');
        set(StatusHandle,'pointer','watch');
        %         set(StatusHandle,'position',statusPos);
        set(StatusHandle,'Menu','none');
        axPos = [.45 .2 .05 .6];
        StatusAxes = axes('Position', axPos, ...
            'XTick',[],...
            'Xlim',[0 1],...
            'Ylim',[0 endNumber],...
            'Parent',StatusHandle,...
            'YColor',FONT_COLOR);
        label = get(StatusAxes,'Xlabel');
        set(label,'string',xlabel);
        label = get(StatusAxes,'Ylabel');
        set(label,'string',ylabel);
        label = get(StatusAxes,'Title');
        set(label,'string',title);
        perCompLine= line('Xdata',[.5 .5],'Ydata',[0 0],...
            'LineWidth',8, 'Color', [1 0 0],...
            'Parent',StatusAxes);
        
        statusData.perCompLine = perCompLine;
        statusData.perComp = 0;
        statusData.endNumber = endNumber;
        set(StatusHandle,'userdata',statusData);
        drawnow;
        
    case 'increment'
        StatusHandle = findobj('Tag','StatusBar');
        statusData = get(StatusHandle(1),'userdata');
        perCompLine = statusData.perCompLine;
        perComp = statusData.perComp;
        if(perComp + incrementAmount < statusData.endNumber)
            perComp = perComp + incrementAmount;
        end
        
        if perComp + incrementAmount == statusData.endNumber
            perComp = statusData.endNumber;
        end
        set(perCompLine,'Ydata',[0 perComp]);
        
        statusData.perComp = perComp;
        set(StatusHandle,'userdata',statusData);
        drawnow;
        
    case 'finished'
        StatusHandle = findobj('Tag','StatusBar');
        Visible = get(StatusHandle, 'visible');
        %         if strcmp('off', Visible)
        %             set(StatusHandle, 'Visible', 'on');
        %             drawnow;
        %         end
        statusData = get(StatusHandle(1),'userdata');
        perCompLine = statusData.perCompLine;
        
        perComp = statusData.endNumber;
        set(perCompLine,'Ydata',[0 perComp]);
        delete(StatusHandle);
        drawnow;
end