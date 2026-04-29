classdef icatb_cls_dice
    % Date 4/20/26
    % Natalia Maksymchuk & Cyrus Eierud
    % Code that performs 
    % Example to calcuate and save state guided ICA:
    %       oc_dice = icatb_cls_dice(dfncInfo);
    %       oc_dice.plot_dice(s_grp1name, s_grp2name, [1,2,5,6,10],
    %              [3,4,7,8,9]); % where subj indices [1,2,5,6,10] may be
    %              a disorder (SZ) and subj indices [3,4,7,8,9] may be control
    %       n_ret = oc_dice.calc_save;

    properties
        dfncInfo
        subj_discrim
    end
    
    methods
        
        function obj = icatb_cls_dice(initialValue)
            % Constructor
            if nargin > 0
                obj.dfncInfo = initialValue;
            end
        end

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

        function G = plot_dice(obj, s_grp0in, s_grp1in, con_subs_grp0, con_subs_grp1)

            % Number of time points / windows, number of correlations, total subjects
            nCorr = length(obj.dfncInfo.comps)^2/2-length(obj.dfncInfo.comps)/2;
            Nsubjects=obj.dfncInfo.userInput.numOfSub;
            Nnetworks = length(obj.dfncInfo.comps);
            
            groupNames.s_grp0 = char(s_grp0in); %make sure they are chars
            groupNames.s_grp1 = char(s_grp1in); %make sure they are chars
               
            for i = 1:Nsubjects
                filename = [obj.dfncInfo.prefix '_sub_' sprintf('%03d', i) '_sess_001_results.mat'];
                temp = load(filename);
                if i ==1
                    N_windows=size(temp.FNCdyn,1);
                    input_dfnc = zeros(N_windows, nCorr, Nsubjects); % Initiate
                end
                input_dfnc(:,:,i) = temp.FNCdyn;
            end
            
            input_dfnc=input_dfnc - min(input_dfnc(:)); %creating tensor with all elements>0
            
            % Computing Dynamic ICE
            for s=1:Nsubjects 

                subConn=diagN(squareform_array(squeeze(input_dfnc(:,:,s))',1));
                for w=1:N_windows 
                    for c=1:Nnetworks
                        h=squeeze(subConn(c,:,w))'; 
                        connDistns(c,:,w,s)=h./nansum(h); % this is the summed connectivity of c to each other network
                        % rescaled to be a distribution.
                        temp_x = log(connDistns(c,:,w,s));
                        temp_z2n = double(temp_x);
                        temp_z2n(find(temp_x==0))=NaN;
                        entropyConnDistns(c,w,s)=-nansum(connDistns(c,:,w,s).*temp_z2n);%entropy of each of these distributions
                        if isreal(entropyConnDistns(c,w,s))==0
                            [c,w,s]
                        end
                    end
                end
            end
            
            meanDynConnEntropy=squeeze(mean(entropyConnDistns,2));%temporal average of entropy of connectivities distributions
            grp0indexes = con_subs_grp0';
            Cindexes = con_subs_grp1';
            Entropygrp0_mean = mean(meanDynConnEntropy(:,grp0indexes), 2); % mean over time and subjects for grp0
            Entropygrp1_mean = mean(meanDynConnEntropy(:,Cindexes), 2); 
            
            Entropygrp0_std = std(meanDynConnEntropy(:,grp0indexes), 0, 2);
            Entropygrp1_std = std(meanDynConnEntropy(:,Cindexes), 0, 2);  % The flag is 0 (default) or 1 to specify normalization by n – 1 or n, 
            % respectively, where n is the number of remaining observations after 
            
            %STD of network connectivity entropies computed across different windows(needed for statistical analysis) 
            STD_DynConnEntropy = squeeze(std(entropyConnDistns,0,2)); % std only across windows - this will be used for statistical analysis
            grp0std_entropy=STD_DynConnEntropy(:,grp0indexes);
            Cstd_entropy=STD_DynConnEntropy(:,Cindexes);
                     
            %% Regression analysis for mean DICE
            B=zeros(Nnetworks,4);
            P=zeros(Nnetworks,4);
            
            alpha=0.05;
            maxflag=1;
            
            demo_fbirn=zeros(length([grp0indexes' Cindexes']),1);
            demo_fbirn(length([grp0indexes])+1:end)=1;

            for c=1:Nnetworks
                 [b1,st1]=robustfit(demo_fbirn, [meanDynConnEntropy(c,grp0indexes) meanDynConnEntropy(c,Cindexes)],'ols');
                 B1(c,:)=b1(2); %for disorder
                 p_u(c,:)=st1.p(2);
                 T1(c,:) = st1.t(2);
            end
            
            p_bhfdr_thresh=fdr(p_u(:),alpha,maxflag);
            
            index=find(p_u(:)<=p_bhfdr_thresh); % FDR correction
            display([num2str(length(index)) ' Networks have p<=' num2str(p_bhfdr_thresh) ' (FDR corrected) :']);
            disp('ICE indexes:')
            obj.dfncInfo.comps(index)
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Tsz=T1(:);
            % Create a mask of zeros
            Tsz_corr = zeros(size(Tsz));
            % Copy values from original data at specified indices
            Tsz_corr(index) = Tsz(index);
            save('DICE_Tval_corr.mat',"Tsz_corr");
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %% Figure for the mean with shaded standard deviation for ICE
            % Define x-axis values (indices)
            x = 1:length(Entropygrp0_mean);
            % Create a figure and plot the mean with shaded standard deviation
            figure('Units', 'normalized', 'Position', [0.1 0.1 0.8 0.4]); % Make figure wider
            hold on;
            % Plot shaded area for grp0 group (in red)
            fill([x, fliplr(x)], [Entropygrp0_mean + Entropygrp0_std; flipud(Entropygrp0_mean - Entropygrp0_std)]', ...
                'r', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', sprintf('SD %s', groupNames.s_grp0));
            % Plot mean line for grp0 group with closed circle markers (in red)
            plot(x, Entropygrp0_mean, '-or', 'LineWidth', 1.5, 'MarkerFaceColor', 'r', 'DisplayName', sprintf('Mean %s', groupNames.s_grp0));
            
            % Plot shaded area for grp1 group (in blue)
            fill([x, fliplr(x)], [Entropygrp1_mean + Entropygrp1_std; flipud(Entropygrp1_mean - Entropygrp1_std)]', ...
                'b', 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'DisplayName', sprintf('SD %s', groupNames.s_grp1));
            % Plot mean line for grp1 group with closed circle markers (in blue)
            plot(x, Entropygrp1_mean, '-ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'DisplayName', sprintf('Mean %s', groupNames.s_grp1));          
            
            % Customize the plot
            xlabel('Networks', 'FontSize', 16, 'FontWeight', 'bold'); % Change x-axis label to 'ICN', make it larger and bold
            ylabel('Mean DICE', 'FontSize', 16, 'FontWeight', 'bold'); % Change y-axis label to 'Mean ICE', make it larger and bold
            if numel(index) > 0
                title(['Mean and SD (shaded) of DICE for ' groupNames.s_grp0 ' and ' groupNames.s_grp1 ', * denotes FDR (\alpha=0.05)'], 'FontSize', 16, 'FontWeight', 'bold');
            else
                title(['Mean and SD (shaded) of DICE for ' groupNames.s_grp0 ' and ' groupNames.s_grp1], 'FontSize', 16, 'FontWeight', 'bold'); % Update title                              
            end

            set(gca, 'XTick', x, 'XTickLabel', obj.dfncInfo.comps)

            % Customize legend (remove box around legend)
            legend({['SD ' groupNames.s_grp0], ['Mean ' groupNames.s_grp0], ['SD ' groupNames.s_grp1], ['Mean ' groupNames.s_grp1]}, 'Location', 'best', 'Box', 'off');            

            % Remove box and grid from figure
            ax = gca;
            ax.Box = 'off';
            ax.XGrid = 'on';
            ax.YGrid = 'on';
            
            % Add star marks for statistically significant networks
            for i = 1:length(index)
                sig_idx = index(i);
                % Get the higher mean between the two groups at this ICN
                max_val = max(Entropygrp1_mean(sig_idx), Entropygrp0_mean(sig_idx));
                % Plot a star a little above the higher mean value
                text(sig_idx, max_val+0.0015, '*', 'FontSize', 30, 'HorizontalAlignment', 'center', 'Color', 'k');
            end

            % Save statistics
            stats_statebased_ttest2.grp0_mean = Entropygrp0_mean;
            stats_statebased_ttest2.grp0_std = Entropygrp0_std;
            stats_statebased_ttest2.grp0_N = length(con_subs_grp0);
            stats_statebased_ttest2.grp1_mean = Entropygrp1_mean;
            stats_statebased_ttest2.grp1_std = Entropygrp1_std;
            stats_statebased_ttest2.grp1_N = length(con_subs_grp1);
            stats_statebased_ttest2.p_u = p_u;
            stats_statebased_ttest2.p_bhfdr_thresh = p_bhfdr_thresh;
            stats_statebased_ttest2.subject_indices.ix_grp0=con_subs_grp0;
            stats_statebased_ttest2.subject_indices.ix_grp1=con_subs_grp1;
            stats_statebased_ttest2.groupNames=groupNames;
            s_stats = fullfile(obj.dfncInfo.outputDir, [obj.dfncInfo.prefix '_stats_dice_results.mat']);
            disp(['DICE summary statistics saved to ' s_stats]);
            icatb_save(s_stats, 'stats_statebased_ttest2');

        end  

        function obj = set_s_outputDir(obj, newValue)
            % Method to set Value
            obj.s_outputDir = newValue;
        end
        
        function val = get_s_outputDir(obj)
            % Method to get Value
            val = obj.s_outputDir;
        end        
    end
end

function Y = squareform_array(X,flag)
    
    if nargin<2
        flag=0;
    end
    
    szX=size(X);
    
    ndims=length(szX);
    
    if flag==0
    if szX(1)~=szX(2)
        error('X must be an array of square matrices with zeros along the diagonal:  first two dimensions should be the same')
    end
    
    if ndims>5
        error('Currently only 3 array dimensions beyond the first two for the square matrix are supported')
    end
    
    if ndims<5
        szX(ndims+1:5)=1;
    end
    
    X=reshape(X,szX);
    
    for d1=1:szX(3)
        for d2=1:szX(4)
            for d3=1:szX(5)
                Y(:,d1,d2,d3)=squareform(diagZ(squeeze(X(:,:,d1,d2,d3))));
            end
        end
    end
    
    Y=sq(Y);
    
    elseif flag==1
    
    if ndims>4
        error('Currently only 3 array dimensions beyond the first one for the vectorized matrix are supported')
    end
    
    if ndims<4
        szX(ndims+1:4)=1;
    end
     
    X=reshape(X,szX);
    
    for d1=1:szX(2)
        for d2=1:szX(3)
            for d3=1:szX(4)
                Y(:,:,d1,d2,d3)=squareform(squeeze(X(:,d1,d2,d3)));
            end
        end
    end
    
    Y=squeeze(Y);
    end
end

function Z = diagN(X)
    %argument X is a scalar, say n.  Z is the nxn matrix with NaN's
    %along the diagonal and zeros everywhere else
    
    szX=size(X);
    
    if length(szX)>5 %| szX(1)~=szX(2)
        error('Must be square matrix or nD n<=5 array of square matrices')
    end
    
    Z=X;
    
    if szX(1)==szX(2)
        if length(szX)<3
            for d=1:szX(1)
                Z(d,d)=NaN;
            end
        elseif length(szX)==3
            for i=1:szX(3)
                for d=1:szX(1)
                    Z(d,d,i)=NaN;
                end
            end
        elseif length(szX)==4
            for i=1:szX(3)
                for j=1:szX(4)
                    
                    for d=1:szX(1)
                        Z(d,d,i,j)=NaN;
                    end
                end
            end
        elseif length(szX)==5
            for i=1:szX(3)
                for j=1:szX(4)
                    for k=1:szX(4)
                        for d=1:szX(1)
                            Z(d,d,i,j,k)=NaN;
                        end
                    end
                end
            end
        end
    end
end

function thresh=fdr(p,alpha,maxflag)
    
    if nargin<2
        alpha=0.05;
    end
    if nargin<3
        maxflag=1;
    end
    
    if mod(sqrt(length(p(:))),1)==0
        p=reshape(p,[sqrt(length(p(:))),sqrt(length(p(:)))]);
    end
    
    if size(p,1)==size(p,2)
        if sum(sum(p-transpose(p)))==0 & (length(unique(diag(p)))==1 | sum(isnan(diag(p)))==size(p,1))%in case of square symmetric matrices with constant diagonals 
            p=squareform(tril(p,-1));
        elseif sum(sum(p-transpose(p)))==0 & (length(unique(diag(p)))>1 | sum(isnan(diag(p)))~=size(p,1))
            p=[squareform(tril(p,-1))';diag(p)];
        end
    end
    
    [pID,pN]=icatb_fdr(p,alpha);
    
    bonthresh=alpha/length(p);
    
    if maxflag==1
        thresh=max([pID,pN,bonthresh]);
    else
        thresh=pID;
    end
end
