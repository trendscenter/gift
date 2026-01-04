classdef icatb_cls_sgica
    % Date 9/2/25
    % Cyrus Eierud
    % Code that performs 
    % Example to calcuate and save state guided ICA:
    %       oc_sgica = icatb_cls_sgica(dfncInfo);
    %       oc_sgica = oc_sgica.set_s_outputDir(outputDir);
    %       oc_sgica = oc_sgica.set_input_FNCdynflat(FNCdynflat);
    %       n_ret = oc_sgica.calc_save;
    % Example to return the report data in structure variable:
    %       oc_sgica = icatb_cls_sgica(dfncInfo);
    %       oc_sgica = oc_sgica.set_s_outputDir(outputDir); 
    properties
        dfncInfo
        input_FNCdynflat
        s_outputDir
    end
    
    methods
        function n_ret = calc_save(obj)
           % try
            % Perform ICA reference guided spatial dfnc
            ncomps = obj.dfncInfo.postprocess.tag_edt_stateguided_numcomps; % ica components

            FNCdynflat_tmp = cat(3, obj.input_FNCdynflat{:});
            FNCdynflat_tmp = permute(FNCdynflat_tmp, [3, 1, 2]);
            Nwin = size(obj.input_FNCdynflat{1}, 1);
            FNCdynflat_tmp = reshape(FNCdynflat_tmp, length(obj.dfncInfo.outputFiles)*Nwin, size(FNCdynflat_tmp, 3));
        
            [whitesig, dewhiteM] = icatb_calculate_pca(FNCdynflat_tmp', ncomps);
            clear FNCdynflat_tmp;
            [~, W_blindICA, A_blindICA, sources] = icatb_icaAlgorithm(1, whitesig'); %gettig template (blindly)
            A_blindICA = dewhiteM*A_blindICA;
            W_blindICA = pinv(A_blindICA);
            priors = sources';
        
            sgica.br_all_sub_fnc_A = zeros(size([obj.input_FNCdynflat{1}; obj.input_FNCdynflat{2}],1)/2,ncomps, obj.dfncInfo.userInput.numOfSess, obj.dfncInfo.userInput.numOfSub);
            sgica.br_all_sub_fnc_S = zeros(ncomps,size(priors,1), obj.dfncInfo.userInput.numOfSess, obj.dfncInfo.userInput.numOfSub);
            % Back Reconstruction constrained_source
            for subject = 1:obj.dfncInfo.userInput.numOfSub
                for nSess = 1:obj.dfncInfo.userInput.numOfSess
                    % Get the matrix from the current cell
                    dfnc = obj.input_FNCdynflat{(subject-1)*obj.dfncInfo.userInput.numOfSess+nSess}; 
                
                    % MOO-ICAR (constrained ICA), A = loadings and sources = components
                    [~, W, sgica.br_all_sub_fnc_A(:,:,nSess,subject), sgica.br_all_sub_fnc_S(:,:,nSess,subject)] = icatb_icaAlgorithm('MOO-ICAR', dfnc, {'ref_data', priors'});
                end
            end
            post_process_file = fullfile(obj.s_outputDir,  [obj.dfncInfo.prefix, '_post_process.mat']);
            sgica.priors=priors;
            save(post_process_file, '-nocompression', 'sgica');

            % Calibration
    
            % Init
            sgica.calib_all_sub_fnc_A = zeros(size(sgica.br_all_sub_fnc_A));
            sgica.calib_all_sub_fnc_S = zeros(size(sgica.br_all_sub_fnc_S));   
            sgica.beta_subj_sess = zeros(obj.dfncInfo.postprocess.tag_edt_stateguided_numcomps, obj.dfncInfo.userInput.numOfSess*obj.dfncInfo.userInput.numOfSub);
    
            num_individuals = (obj.dfncInfo.userInput.numOfSub-1)*obj.dfncInfo.userInput.numOfSess+obj.dfncInfo.userInput.numOfSess;
            % Process each individual within the group
            for subject = 1:obj.dfncInfo.userInput.numOfSub
                for nSess = 1:obj.dfncInfo.userInput.numOfSess    
                    ix_scan = (subject-1)*obj.dfncInfo.userInput.numOfSess+nSess;
                    % Load dFNC data for the current individual
                    dFNC_ind = obj.input_FNCdynflat{(subject-1)*obj.dfncInfo.userInput.numOfSess+nSess}; % 
                    % Load constrained sourcesfor the current scan
                    predictors = zeros(numel(dFNC_ind), obj.dfncInfo.postprocess.tag_edt_stateguided_numcomps);
                    % Prepare the predictors by computing a_i * s_i for each component i
                    for i1 = 1:obj.dfncInfo.postprocess.tag_edt_stateguided_numcomps
                        % Extract component a_i (windows-by-1 vector)
                        a_i = sgica.br_all_sub_fnc_A(:, i1, nSess, subject); % 
                        s_i = sgica.br_all_sub_fnc_S(i1, :, nSess, subject); % 
                        % Compute the product a_i * s_i to create the predictor for component i
                        component_vector = a_i * s_i; % 115 x 1378
                        predictors(:, i1) = component_vector(:); %flatten & matrix
                    end
                    response = dFNC_ind(:); % Flatten dFNC data
%remove mean or zscore ?

                    beta = regress(response, predictors);
                    % Apply beta weights and scale each component by std
                    for i2 = 1:obj.dfncInfo.postprocess.tag_edt_stateguided_numcomps
                        % Scale the component (sica_calib.fnc_S)
                        std_i = std(sgica.br_all_sub_fnc_S(i2, :, nSess, subject)); % Compute standard deviation of the scaled component
                        sgica.calib_all_sub_fnc_S(i2, :, nSess, subject) = sgica.br_all_sub_fnc_S(i2, :, nSess, subject)/std_i;
                        % sica_calib.fnc_A = final_timecourse + (sica_calib.fnc_S / std_i) * std(A_ind(:, i));
                        sgica.calib_all_sub_fnc_A(:,i2, nSess, subject) = beta(i2) * sgica.br_all_sub_fnc_A(:, i2, nSess, subject) * std_i;               
                    end
                end
            end

            % Create z-scored calibrated results
            sgica.calib_all_sub_fnc_A_z = zeros(size(sgica.calib_all_sub_fnc_A));
            for n_ic=1:size(sgica.calib_all_sub_fnc_A_z,2)
                for n_sess=1:size(sgica.calib_all_sub_fnc_A_z,3)
                    for n_subj=1:size(sgica.calib_all_sub_fnc_A_z,4)
                        % Just Time courses for now
                        sgica.calib_all_sub_fnc_A_z(:,n_ic, n_sess, n_subj) = icatb_zscore(sgica.calib_all_sub_fnc_A(:,n_ic, n_sess, n_subj));
                    end
                end
            end

            save(post_process_file, '-nocompression', '-append', 'sgica');
            n_ret=0; % everything ran through
     %   catch
      %      n_ret=1; % something failed
     %   end
        end

        function con_dwell_all_flat = ar_report_dwell(obj, con_dwell_all)
            con_dwell_all_flat = permute(con_dwell_all,[3,1,2]);
            con_dwell_all_flat = reshape(con_dwell_all_flat,[], size(con_dwell_all_flat,2)*size(con_dwell_all_flat,3));
            con_dwell_all_flat = con_dwell_all_flat';
            %calib_all_sub_fnc_A_absmax = max(abs(calib_all_sub_fnc_A_flat));
        end

        function G = plot_mn_dwell(obj,con_dwell_all_flat, G, n_windows)
            % function G = icatb_dfnc_plot_statevector_stats(k, F, TM, MDT, NT, G)
            k=size(con_dwell_all_flat,2);

            if (~exist('G', 'var'))
                G = figure;
            end            
            
            if isvector(con_dwell_all_flat)
                SS = 1; % single subject
            else
                SS = 0; % group
            end
            
            subplot(2,1,1)
            if SS
                plot(1:k, con_dwell_all_flat./n_windows, 'm');
            else
                icatb_plot_with_ste_area(gca, 1:k, con_dwell_all_flat./n_windows, [], 'm');
            end
            box off; set(gca, 'TickDir', 'out')
            ylabel('Frequency')
            xlabel('Reference Guided Spatial dFNC (ICA index)')
            
            
            subplot(2,1,2)
            if SS
                plot(1:k, con_dwell_all_flat, 'm');
            else
                icatb_plot_with_ste_area(gca, 1:k, con_dwell_all_flat, [], 'm');
            end
            box off; set(gca, 'TickDir', 'out')
            ylabel('Mean dwell time (windows)')
            xlabel('Reference Guided Spatial dFNC (ICA index)')
            
        end

        function G = plot_transitions(obj,con_tran_pos_flat,con_tran_neg_flat, G)
            k=size(con_tran_pos_flat,2);

            if (~exist('G', 'var'))
                G = figure;
            end            
            
            if isvector(con_tran_pos_flat)
                SS = 1; % single subject
            else
                SS = 0; % group
            end
            
            subplot(2,1,1)
            if SS
                plot(1:k, con_tran_pos_flat, 'm');
            else
                icatb_plot_with_ste_area(gca, 1:k, con_tran_pos_flat, [], 'm');
            end
            box off; set(gca, 'TickDir', 'out')
            ylabel('Transitions into Positive percentile')
            xlabel('Reference Guided Spatial dFNC (ICA index)')
            
            
            subplot(2,1,2)
            if SS
                plot(1:k, con_tran_neh_flat, 'm');
            else
                icatb_plot_with_ste_area(gca, 1:k, con_tran_neg_flat, [], 'm');
            end
            box off; set(gca, 'TickDir', 'out')
            ylabel('Transitions into Negative percentile')
            xlabel('Reference Guided Spatial dFNC (ICA index)')
            
        end

        function obj = cls_stat_two_sample(obj, figData)

            disp('Selected design criteria is two sample t-test for Reference Guided Spatial dFNC');
            

            %CE091625 CREATE THE ARRAY (replace below by creating correlations in accordance with Vicne meeting 091525 - fold6im)
            % VARIABLE ar_report_dwell
%             %% Two sample t-test
%             if (size(dfnc_corrs, 1) > 1)
%                 dfnc_corrs = squeeze(mean(dfnc_corrs, 1));
%             else
%                 dfnc_corrs = squeeze(dfnc_corrs);
%             end
            
%        %CE091624 DONT know
%             if (exist('dfncTaskConnectivity', 'var'))
%                 if (size(dfncTaskConnectivity, 1) > 1)
%                     dfncTaskConnectivity = squeeze(mean(dfncTaskConnectivity, 1));
%                 else
%                     dfncTaskConnectivity = squeeze(dfncTaskConnectivity);
%                 end
%             end
            
            % numClusters = size(dfnc_corrs, 3);
% % % %  ce121725            numClusters = size(ar_report_dwell, ?);
            
            %% Initialize results
            t_u = cell(1, numClusters);
            p_u = t_u;
            stats_u = t_u;
            N = zeros(2, numClusters);
            mean_u = cell(2, numClusters);
            
            g1 = [figData.ttestOpts.t.val{1}];
            g2 = [figData.ttestOpts.t.val{2}];
            subject_indices = cell(2, numClusters);
            
%             %% Compute and save
%             for nC = 1:numClusters
%                 disp(['Computing two sample t-test on cluster state# ', num2str(nC), ' ...']);
%                 %tmp = squeeze(dfnc_corrs(:, :, nC));
%                 tmp1 =  squeeze(dfnc_corrs(g1, :, nC));
%                 tmp2 =  squeeze(dfnc_corrs(g2, :, nC));
%                 
%                 chk1 = find(isfinite(tmp1(:, 1)) == 1);
%                 chk2 = find(isfinite(tmp2(:, 1)) == 1);
%                 
%                 if (~isempty(chk1) && ~isempty(chk2))
%                     if ((length(chk1) + length(chk2)) > 2)
%                         tmp1 = tmp1(chk1, :);
%                         tmp2 = tmp2(chk2, :);
%                         N1 = length(chk1);
%                         N2 = length(chk2);
%                         modelX = ones(N1 + N2, 1);
%                         modelX(N1 + 1:end) = 0;
%                         N(1, nC) = N1;
%                         N(2, nC) = N2;
%                         subject_indices{1, nC} = g1(chk1);
%                         subject_indices{2, nC} = g2(chk2);
%                         mean_u{1, nC} = mean(tmp1);
%                         mean_u{2, nC} = mean(tmp2);
%                         [t_u{nC}, p_u{nC}, stats_u{nC}] = mT([tmp1;tmp2], modelX, [], 1, {'verbose'});
%                     end
%                 end
%                 
%             end
%             
%             outFile = fullfile(cluster_stats_directory, [figData.prefix, '_two_sample_ttest_results.mat']);
%             disp(['Two sample t-test results are saved in ', outFile]);
%             fprintf('\n\n');
%             save(outFile, 't_u', 'p_u', 'stats_u', 'mean_u', 'N', 'groupNames', 'groupVals', 'subject_indices');
%             
%             
%             if (exist('dfncTaskConnectivity', 'var'))
%                 
%                 
%                 modelX = ones(length(g1) + length(g2), 1);
%                 modelX(length(g1) + 1:end) = 0;
%                 
%                 clear t_task_conn p_task_conn stats_task_conn mean_task_conn
%                 
%                 mean_task_conn = cell(2, length(selectedRegressors));
%                 t_task_conn = cell(1, length(selectedRegressors));
%                 p_task_conn = cell(1, length(selectedRegressors));
%                 stats_task_conn = cell(1, length(selectedRegressors));
%                 
%                 for nRegress = 1:length(selectedRegressors)
%                     disp(['Computing one sample t-test on task connectivity (', selectedRegressors{nRegress}, ' ...']);
%                     tmp1 = squeeze(dfncTaskConnectivity(g1, :, nRegress));
%                     tmp2 = squeeze(dfncTaskConnectivity(g2, :, nRegress));
%                     mean_task_conn{1, nRegress} = mean(tmp1);
%                     mean_task_conn{2, nRegress} = mean(tmp2);
%                     [t_task_conn{nRegress}, p_task_conn{nRegress}, stats_task_conn{nRegress}] = mT([tmp1;tmp2], modelX, [], 1, {'verbose'});
%                 end
%                 save(outFile, 't_task_conn', 'p_task_conn', 'stats_task_conn', 'mean_task_conn', '-append');
%             end

        end

        
        function obj = icatb_cls_sgica(initialValue)
            % Constructor
            if nargin > 0
                obj.dfncInfo = initialValue;
            end
        end
        
        function obj = set_dfncInfo(obj, newValue)
            % Method to set Value
            obj.dfncInfo = newValue;
        end
        
        function val = get_dfncInfo(obj)
            % Method to get Value
            val = obj.dfncInfo;
        end

        function obj = set_s_outputDir(obj, newValue)
            % Method to set Value
            obj.s_outputDir = newValue;
        end
        
        function val = get_s_outputDir(obj)
            % Method to get Value
            val = obj.s_outputDir;
        end        

        function obj = set_input_FNCdynflat(obj, newValue)
            % Method to set Value
            obj.input_FNCdynflat = newValue;
        end

        function val = get_input_FNCdynflat(obj)
            % Method to get Value
            val = obj.input_FNCdynflat;
        end
    end
end
