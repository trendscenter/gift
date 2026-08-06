function mat_colormap = trd_fnc_cm_jetwhite(Bar_range, n_colors)
        % Method to set Value
%             Bar_range=[-1,1]; %use the bar range [-1,1] for correlation
%             
%             % TReND colormap
%             % Bar_range is two element vector of lowest and highest val in colormap
%             % example: mymap=fun_trends_cmap=[-1,1] %for correlation
%             % blue to white middle and then to red colormap

        % Default: do not reduce the colormap
        if nargin < 2
            b_n_colors = false;
        end

        %% gradient map
        ini_idx = 20; % small:edge close to black;
        zero_modified_idx = 2;
        ratio = (Bar_range(2) - 0)/(0 - Bar_range(1));
        num_hot = 100;
        num_col = round((num_hot-zero_modified_idx)/ratio);
        % hot map
        colormap(rand(num_hot,3))
        hot_map = colormap(hot);
        hot_map_use = hot_map(end-zero_modified_idx:-1:(1+ini_idx),:);
        % col map
        colormap(rand(num_col,3))
        hot_map = colormap(hot);
        col_map_use = zeros(size(hot_map,1)-round(ini_idx/ratio),size(hot_map,2));
        col_map_use(:,1) = hot_map((1+round(ini_idx/ratio):end),3);
        col_map_use(:,2) = hot_map((1+round(ini_idx/ratio):end),2);
        col_map_use(:,3) = hot_map((1+round(ini_idx/ratio):end),1);
        
        mat_colormap = [col_map_use;hot_map_use];

        % Reduce to n_colors if requested
        if b_n_colors
            xOld = linspace(0,1,size(mat_colormap,1));
            xNew = linspace(0,1,n_colors);
            mat_colormap = interp1(xOld, mat_colormap, xNew, 'linear');
        end        
end