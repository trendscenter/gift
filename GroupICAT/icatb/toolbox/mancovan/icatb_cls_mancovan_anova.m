classdef icatb_cls_mancovan_anova
    % Date 10/27/25
    % Cyrus Eierud
    % Code that implements ANOVA into MANCOVAN
    % Example to return the stats type:
    %       oc_mancovan_anova = icatb_cls_mancovan_anova('1-way, x-level anova');
    %       desCriteria = oc_mancovan_anova.get_s_desCriteria();     
    properties
        s_desCriteria
    end
    
    methods
        function oc = icatb_cls_mancovan_anova(initialValue)
            % Constructor
            if nargin > 0
                oc.s_desCriteria = initialValue;
            end
        end

        function val = get_s_desCriteria(oc)
            % Method to get Value
            val = oc.s_desCriteria;
        end

        function [cs_level_names ccoi_level_subs] = get_mcs_levels_GUI(oc, subjectStr)
            % Method to get anova levels
            n_levels = icatb_inputdlg2('Enter number of ANOVA levels', 'ANOVA Levels', 1, {''});
            n_levels = round(str2num(n_levels{1}));
            cs_level_names = cell(n_levels,1);
            ccoi_level_subs = cell(n_levels,1); 
        
            for n_level = 1:n_levels
                [cs_level_names{n_level}, ccoi_level_subs{n_level}] = icatb_select_groups_gui(subjectStr, ['Group ' num2str(n_level)], ['selGrp  ' num2str(n_level)], '', []);
            end
        end

        function [stats_u, t_u, p_u, tmp_con_name] = m_calc_anova(oc, data, mancovanInfo)
            [i_subs, i_features] = size(data);
            
            coi_group_subs = zeros(i_subs,1);
            for i_group = 1:length(mancovanInfo.userInput.ttestOpts.anova.val{1, 1})
                ix_grp = mancovanInfo.userInput.ttestOpts.anova.val{1, 1}{i_group,1};
                coi_group_subs(ix_grp) = i_group;
            end
            
            % Preallocate outputs (store F in t_u{1} to mimic your interface)
            F_u{1}     = NaN(1, i_features);   % will hold F-statistic for each column
            p_u{1}     = NaN(1, i_features);   % p-value per column
            stats_u{1}.Terms = cell(1,1);
            %stats_u{1}.X = NaN(1, i_features);
            stats_u{1}.SSE = NaN(1, i_features);
            stats_u{1}.DFE = NaN(1,i_features);
            stats_u{1}.MSE = NaN(1, i_features);
            %stats_u{1}.B = NaN(1, i_features);
            stats_u{1}.Levels = [1,0];
            stats_u{1}.Term = 1;           % scalar => treated as "main effect"
            
            for j = 1:i_features
                y = data(:, j);
                % handle NaNs (optional; anova1 ignores NaNs in y, but we match group length)
                idx = ~isnan(y);
                if sum(idx) < 3 || numel(unique(coi_group_subs(idx))) < 2
                    continue  % not enough data/levels
                end
            
                [p, tbl, stats_an] = anova1(y(idx), coi_group_subs(idx), 'off');  % no figure
            
                % Extract ANOVA pieces
                Fval = tbl{2,5};      % F
                SS   = tbl{2,2};      % treatment sum of squares
                MS   = tbl{2,4};      % treatment mean square
            
                % Fill "like" your [t_u, p_u, stats_u] interface
                t_u{1}(j)          = Fval;         % using F in place of t
                p_u{1}(j)          = p;
                stats_u{1}.F(j)    = Fval;
                stats_u{1}.DFE(j)   = stats_an.df; 
                stats_u{1}.SSE(j)   = SS;
                stats_u{1}.MSE(j)   = MS;
                
    %             if ~isfield(stats_u{1}(j),'Levels') || isempty(stats_u{1}(j).Levels)
    %                 stats_u{1}(j).Levels = [1 0]; % minimal placeholder so get_contrast_label doesn't error
    %             end
            end
            if length(mancovanInfo.userInput.ttestOpts.anova.val{1, 1}) == 2
                tmp_con_name = ['ANOVA Levels: ' mancovanInfo.userInput.ttestOpts.anova.name{1, 1}{1,1} ', ' mancovanInfo.userInput.ttestOpts.anova.name{1, 1}{2,1}];
            elseif length(mancovanInfo.userInput.ttestOpts.anova.val{1, 1}) == 3
                tmp_con_name = ['ANOVA Levels: ' mancovanInfo.userInput.ttestOpts.anova.name{1, 1}{1,1} ', ' mancovanInfo.userInput.ttestOpts.anova.name{1, 1}{2,1} ', ' mancovanInfo.userInput.ttestOpts.anova.name{1, 1}{3,1}];
            else
                tmp_con_name = [num2str(length(mancovanInfo.userInput.ttestOpts.anova.val{1, 1})) ' ANOVA Levels'];
            end            
        end

        function [stats_u, F_u, p_u, tmp_con_name, N, subject_indices, mean_u] = m_calc_anova_dfnc(oc, data, groupVals, groupNames)
            [i_subs, i_features, i_states] = size(data);
            
            coi_group_subs = zeros(i_subs,1);
            for i_group = 1:length(groupNames)
                ix_grp = groupVals{i_group,1};
                coi_group_subs(ix_grp) = i_group;
            end


            % Preallocate outputs (store F in t_u{1} to mimic your interface)
% %             %F_u{1}     = NaN(i_features, i_states);   % will hold F-statistic for each column
            F_u = cell(1, i_states);
% %             % p_u{1}     = NaN(i_features, i_states);   % p-value per column
            p_u = cell(1, i_states);
% %             %stats_u{1}.Terms = cell(1, i_states);
            stats_u = cell(1, i_states);
            %stats_u{1}.X = NaN(1, i_features);
% % % %             stats_u{1}.SSE = NaN(i_features, i_states);
% % % %             stats_u{1}.DFE = NaN(i_features, i_states);
% % % %             stats_u{1}.MSE = NaN(i_features, i_states);
            %stats_u{1}.B = NaN(1, i_features);
            %stats_u{1}.Levels = [1,0];
            %stats_u{1}.Term = 1;           % scalar => treated as "main effect"

            for i_state = 1:i_states     
                disp(['Computing ANOVA on cluster state# ', num2str(i_state), ' ...']);

                ar_tmp_grp = zeros(length(find(coi_group_subs==i_group)),1);
                clear ccod_tmp_grp c_chk_grp coi_tot;
                coi_tot =[];
                for i_group = 1:length(groupNames)
                    ccod_tmp_grp{i_group,1} = squeeze(data(find(coi_group_subs == i_group), :, i_state));
                    c_chk_grp{i_group,1} = find(isfinite(ccod_tmp_grp{i_group,1}(:, 1)) == 1);              
                    coi_tot=[coi_tot ; c_chk_grp{i_group,1}];
                end

                if (~isempty(coi_tot))
                    if (length(coi_tot) > 2)
                        for i_group = 1:length(groupNames)
                            ccod_tmp_grp{i_group,1} = ccod_tmp_grp{i_group,1}(c_chk_grp{i_group,1}, :);
                            N_tmp{i_group,1} = size(ccod_tmp_grp{i_group,1},1);
                            N{i_group,i_state} = N_tmp{i_group,1};
                            ix_current = find(coi_group_subs == i_group);
                            subject_indices{i_group, i_state} = ix_current(c_chk_grp{i_group,1});
                            mean_u{i_group, i_state} = mean(ccod_tmp_grp{i_group,1});
                            %ce102925 [t_u{nC}, p_u{nC}, stats_u{nC}] = mT([tmp1;tmp2], modelX, [], 1, {'verbose'});
                        end
%                         modelX = ones(N1 + N2, 1);
%                         modelX(N1 + 1:end) = 0;
                    end
                end                   
                F_u_tmp=NaN(1, i_features);
                p_u_tmp=NaN(1, i_features);
                F_tmp=NaN(1, i_features);
                DFE_tmp=NaN(1, i_features);
                SSE_tmp=NaN(1, i_features);
                MSE_tmp=NaN(1, i_features);                  
                for j = 1:i_features
                    y = squeeze(data(:, j, i_state));
                    % handle NaNs (optional; anova1 ignores NaNs in y, but we match group length)
                    idx = ~isnan(y);
                    if sum(idx) < 3 || numel(unique(coi_group_subs(idx))) < 2
                        continue  % not enough data/levels
                    end
                
                    [p, tbl, stats_an] = anova1(y(idx), coi_group_subs(idx), 'off');  % no figure
                
                    % Extract ANOVA pieces
                    Fval = tbl{2,5};      % F
                    SS   = tbl{2,2};      % treatment sum of squares
                    MS   = tbl{2,4};      % treatment mean square
                
% % %                     % Fill "like" your [t_u, p_u, stats_u] interface
% % %                     F_u{1}(j, i_state)          = Fval;         % using F in place of t
% % %                     p_u{1}(j, i_state)          = p;
% % %                     stats_u{1}.F(j, i_state)    = Fval;
% % %                     stats_u{1}.DFE(j, i_state)   = stats_an.df; 
% % %                     stats_u{1}.SSE(j, i_state)   = SS;
% % %                     stats_u{1}.MSE(j, i_state)   = MS;
                    F_u_tmp(1,j) = Fval;% using F in place of t    
                    p_u_tmp(1,j) = p;
                    F_tmp(1,j) = Fval;
                    DFE_tmp(1,j) = stats_an.df;
                    SSE_tmp(1,j) = SS;
                    MSE_tmp(1,j) = MS;
        %             if ~isfield(stats_u{1}(j),'Levels') || isempty(stats_u{1}(j).Levels)
        %                 stats_u{1}(j).Levels = [1 0]; % minimal placeholder so get_contrast_label doesn't error
        %             end
                end
                F_u{i_state} = F_u_tmp;         
                p_u{i_state} = p_u_tmp;
                stats_u{i_state}.F = F_tmp;      
                stats_u{i_state}.DFE = DFE_tmp; 
                stats_u{i_state}.SSE = SSE_tmp;
                stats_u{i_state}.MSE = MSE_tmp;     
                clear Fval_tmp p_u_tmp F_tmp DFE_tmp SSE_tmp MSE_tmp;
            end
            if length(groupVals) == 2
                tmp_con_name = ['ANOVA Levels: ' groupNames{1, 1} ', ' groupNames{2, 1}];
            elseif length(groupVals) == 3
                tmp_con_name = ['ANOVA Levels: ' groupNames{1, 1} ', ' groupNames{2, 1} ', ' groupNames{3, 1}];
            else
                tmp_con_name = [num2str(length(groupVals{1, 1})) ' ANOVA Levels'];
            end

        end        
    end
end

