classdef icatb_dfnc_stats_cls
    % Date 12/31/25
    % Cyrus Eierud
    % Code that performs statistics for icatb_dfnc_stats
    % Example to calcuate two sample ttest on correlations:
    %       oc_sgica = icatb_dfnc_stats_cls();
    %       [t_u, p_u, stats_u, mean_u, N, subject_indices] = oc_sgica.m_ttest2(in_array_subj_vals_states, ix_grp1, ix_grp2);
%     properties
%         s
%     end
    
    methods

        function [t_u, p_u, stats_u, mean_u, N, subject_indices] = m_ttest2(obj, statvals_subjxvalsxclusters, grp1, grp2)
            % cls_m_ttest2 does a two sample ttest of array statvals_subjxvalsxclusters
            % grp1 is the column vector of subject indexes for grp1
            % grp2 is the column vector of subject indexes for grp2
            disp('Design criteria two sample t-test will be performed');

            numClusters = size(statvals_subjxvalsxclusters, 3);
            
            %% Initialize results
            t_u = cell(1, numClusters);
            p_u = t_u;
            stats_u = t_u;
            N = zeros(2, numClusters);
            mean_u = cell(2, numClusters);
    
            subject_indices = cell(2, numClusters);
            
            %% Compute and save
            for nC = 1:numClusters
                disp(['Computing two sample t-test on cluster state# ', num2str(nC), ' ...']);
                tmp1 =  squeeze(statvals_subjxvalsxclusters(grp1, :, nC));
                tmp2 =  squeeze(statvals_subjxvalsxclusters(grp2, :, nC));
                
                chk1 = find(isfinite(tmp1(:, 1)) == 1);
                chk2 = find(isfinite(tmp2(:, 1)) == 1);
                
                if (~isempty(chk1) && ~isempty(chk2))
                    if ((length(chk1) + length(chk2)) > 2)
                        tmp1 = tmp1(chk1, :);
                        tmp2 = tmp2(chk2, :);
                        N1 = length(chk1);
                        N2 = length(chk2);
                        modelX = ones(N1 + N2, 1);
                        modelX(N1 + 1:end) = 0;
                        N(1, nC) = N1;
                        N(2, nC) = N2;
                        subject_indices{1, nC} = grp1(chk1);
                        subject_indices{2, nC} = grp2(chk2);
                        mean_u{1, nC} = mean(tmp1);
                        mean_u{2, nC} = mean(tmp2);
                        [t_u{nC}, p_u{nC}, stats_u{nC}] = mT([tmp1;tmp2], modelX, [], 1, {'verbose'});
                    end
                end
                
            end


        end

        
        function obj = icatb_cls_sgica()
            % Constructor does nothing
            obj = [];
        end
        
    end
end
