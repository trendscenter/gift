function icatb_dfnc_results_plot_buttons(graphicsH)
    % number of figures
    xOffset = 0.06;
    yOffset = 0.01;
    buttonWidth = 0.18;
    buttonHeight = 0.025;
    xOrigin = 1 - 3*xOffset - 3*buttonWidth;

    numFigures = length(graphicsH);

    for numFig = numFigures:-1:1
        figureData = get(graphicsH(numFig).H, 'userdata');
        figure(graphicsH(numFig).H);
        figureData.index = numFig;
        figureData.numFigures = numFigures;
        figureData.GraphicsHandle = graphicsH;

        % --Buttons
        % back button and position
        BackButtonPos = [xOrigin yOffset buttonWidth buttonHeight];
        BackButton = icatb_uicontrol('parent', graphicsH(numFig).H, 'style', 'pushbutton', 'units', ...
            'normalized', 'position', BackButtonPos, 'string', 'Ultra Contrast', 'callback', ...
            {@FigureCallback_contrastmore, graphicsH(numFig).H});
        NextButtonPos =[BackButtonPos(1) + BackButtonPos(3) + xOffset BackButtonPos(2) buttonWidth buttonHeight];
        NextButton = icatb_uicontrol('parent', graphicsH(numFig).H, 'style', 'pushbutton', 'units', 'normalized', ...
            'position', NextButtonPos, 'string', 'Normal Contrast', 'callback', {@FigureCallback_contrast_norm, graphicsH(numFig).H});

        try
            set(graphicsH(numFig).H, 'toolbar', 'figure');
        catch
        end
        clear figureData;
    end
end

% previous figure callback
function FigureCallback_contrast_norm(handleObj, event_data, handles)

    h_ax_tmp = handles.CurrentAxes;
    d_color_max = max(h_ax_tmp.Colorbar.Limits);
    d_color_min = min(h_ax_tmp.Colorbar.Limits);
    oc_fnc = icatb_cls_fnc_misc();
    if d_color_min >= 0
        %only upper half
        handles.Colormap = oc_fnc.get_colors_jet_white([0,1]);
    else
        %from middle of colors half
        handles.Colormap = oc_fnc.get_colors_jet_white([-1,1]);
    end
end


% next figure callback
function FigureCallback_contrastmore(handleObj, event_data, handles)

    h_ax_tmp = handles.CurrentAxes;
    d_color_max = max(h_ax_tmp.Colorbar.Limits);
    d_color_min = min(h_ax_tmp.Colorbar.Limits);
    oc_fnc = icatb_cls_fnc_misc();
    if d_color_min >= 0
        %only upper half
        mat_col_jet_white = oc_fnc.get_colors_jet_white([0,1]);
        mat_col_jet_white_xtra = mat_col_jet_white(end,:);
        for i_cols=1:round(length(mat_col_jet_white)*.65)
            mat_col_jet_white_xtra = [mat_col_jet_white_xtra;mat_col_jet_white_xtra];
        end
        handles.Colormap = [mat_col_jet_white;mat_col_jet_white(end,:)];
    else
        mat_col_jet_white = oc_fnc.get_colors_jet_white([-1,1]);
        mat_col_jet_white_cold = mat_col_jet_white(1,:);
        mat_col_jet_white_hot = mat_col_jet_white(end,:);
        for i_cols=1:round(length(mat_col_jet_white)*0.65/2)
            mat_col_jet_white_cold = [mat_col_jet_white_cold;mat_col_jet_white(1,:)];
            mat_col_jet_white_hot = [mat_col_jet_white_hot;mat_col_jet_white(end,:)];
        end
        handles.Colormap = [mat_col_jet_white_cold;mat_col_jet_white;mat_col_jet_white_hot];
    end

end