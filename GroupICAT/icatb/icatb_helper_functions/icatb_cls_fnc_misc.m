classdef icatb_cls_fnc_misc

%%%%%%%%%%%%%DELETE DELETE

    % Date 12/31/25
    % Cyrus Eierud
    % Code that performs statistics for icatb_dfnc_stats
    % Example to calcuate two sample ttest on correlations:
    %       oc_sgica = icatb_cls_fnc_misc();
    %       [t_u, p_u, stats_u, mean_u, N, subject_indices] = oc_sgica.get_colors_jet_white(in_array_subj_vals_states, ix_grp1, ix_grp2);
%     properties
%         s
%     end
    
    methods

        function mat_colormap = get_colors_jet_white(obj, Bar_range)
            % Method to set Value
%             Bar_range=[-1,1]; %use the bar range [-1,1] for correlation
%             
%             % TReND colormap
%             % Bar_range is two element vector of lowest and highest val in colormap
%             % example: mymap=fun_trends_cmap=[-1,1] %for correlation
%             % blue to white middle and then to red colormap
            
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
        end
        
        function obj = icatb_cls_sgica()
            % Constructor does nothing
            obj = [];
        end
        
    end
end
