function varargout = icatb_fb_results_plot(dfncInfo, display_criteria, statsInfo)
    %% DFNC results
    %

    icatb_defaults;
    global FONT_COLOR;
    global UI_FONTNAME;
    global UI_FS;
    global DFNC_DEFAULTS;

    if (~exist('display_criteria', 'var'))
        display_criteria = 'fnc oscillations';
    end

    if (~exist('statsInfo', 'var'))
        statsInfo.threshdesc = 'none';
        statsInfo.p_threshold = 0.05;
        statsInfo.outputDir = pwd;
    end

    comps = [dfncInfo.comps];
    comps = comps(:);

    outputDir = dfncInfo.outputDir;

    cd(outputDir);

    post_process_file = fullfile(outputDir,  [dfncInfo.prefix, '_post_process.mat']);

    load(post_process_file);

    figVisible = 'on';
    if (nargout == 1)
        figVisible = 'off';
        if (~strcmpi(display_criteria, 'group comparisons'))
            html_dir = fullfile(dfncInfo.outputDir, 'html');
            if (~exist(html_dir, 'dir'))
                mkdir(outputDir, 'html');
            end
        else
            cluster_stats_directory = statsInfo.outputDir;
            html_dir = fullfile(statsInfo.outputDir, 'html');
            if (~exist(html_dir, 'dir'))
                mkdir(html_dir);
            end
        end
    end

    network_values = zeros(1, length(dfncInfo.userInput.comp));
    for nV = 1:length(network_values)
        network_values(nV) = length(dfncInfo.userInput.comp(nV).value);
    end
    network_names =  cellstr(char(dfncInfo.userInput.comp.name));

    if (length(network_names) == 1)
        network_names = '';
    end

    comps = [dfncInfo.comps];
    comps = comps(:);

    TR = min(dfncInfo.TR);

    n_win = 1;

    if (strcmpi(display_criteria, 'fnc oscillations'))
        for n_filterbank = 1:size(FNCamp,2)
            
            %% FNC Oscillations
            H(n_win).H = icatb_getGraphics(['FNC Oscillations (1), Band=' num2str(n_filterbank)], 'graphics', 'dfnc_summary', figVisible);
            
            colormap(jet);
            
            %sh = subplot(1, 1, 1);
            sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
            
            CLIM = [min(FNCamp(:,n_filterbank)), max(FNCamp(:,n_filterbank))];
            FNCamp_plot = icatb_vec2mat(FNCamp(:,n_filterbank), 1);
            trd_fnc_plot(FNCamp_plot, CLIM, cellstr(num2str(comps)), (1:length(comps)), H(n_win).H, 'Average Std Of FNC Spectra', ...
                sh, network_values, network_names);
            title(['Average Std Of FNC Spectra Band=' num2str(n_filterbank)], 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            
            if (nargout == 1)
                outFile = [dfncInfo.prefix, '_avg_std_oscillations_band' num2str(n_filterbank) '.jpg'];
                printFile(H(1).H, fullfile(html_dir, outFile));
                results(1).file = outFile;
                results(1).text = ['Spectra is computed on dFNC correlations for band ' num2str(n_filterbank) '. Standard deviation of spectra is then averaged across subjects.'];
            end
            
            n_win = n_win+1;                            
            H(n_win).H = icatb_getGraphics(['FNC Oscillations (2), Band=' num2str(n_filterbank)], 'graphics', 'dfnc_summary', figVisible);
            
            colormap(jet);
            
            sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
            
            CLIM = [min(FNCcm(:,n_filterbank)), max(FNCcm(:,n_filterbank))];
            FNCcm_plot = icatb_vec2mat(FNCcm(:,n_filterbank), 1);
            trd_fnc_plot(FNCcm_plot, CLIM, cellstr(num2str(comps)), (1:length(comps)), H(n_win).H, ['Average Center-Of-Mass Of FNC Spectra Band=' num2str(n_filterbank)], sh, ...
                network_values, network_names);
            title(['Average Center-Of-Mass Of FNC Spectra for band=' num2str(n_filterbank)], 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            
            if (nargout == 1)
                outFile = [dfncInfo.prefix, '_avg_center_of_mass_oscillations_band' num2str(n_filterbank) '.jpg'];
                printFile(H(end).H, fullfile(html_dir, outFile));
                results(2).file = outFile;
                results(2).text = ['Spectral center of mass is computed on dFNC correlations for band=' num2str(n_filterbank) ' and averaged across subjects.'];
            end
            n_win = n_win+1; 
        end
    
    else
        if ~exist('clusterInfo','var')
            error('Cluster stats doesn''t exist. You need to re-run the post-processing step');
        end
        %% Cluster centroids
        %for n_filterbank = 1:size(FNCamp,2)
        M = length (dfncInfo.outputFiles);
        Nwin = length(clusterInfo.IDXall) / M;
        aIND = reshape(clusterInfo.IDXall, M, Nwin);
        aIND = aIND';
        for ii = 1:dfncInfo.postprocess.num_clusters
            fig_title = ['Cluster Centroid (', num2str(ii), ')'];
            H(ii).H = icatb_getGraphics(fig_title, 'graphics', 'dfnc_summary3', figVisible);
            colormap(trd_fnc_cm_jetwhite([-1,1], 64));
            sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
            CLIM = max(abs(clusterInfo.Call(ii,:)));
            CLIM = [-CLIM, CLIM];
            tmp = icatb_vec2mat(clusterInfo.Call(ii, :), 1);
            titleStr = sprintf('%d (%d%%)', sum(clusterInfo.IDXall == ii),  round(100*sum(clusterInfo.IDXall == ii)/length(clusterInfo.IDXall)));
            trd_fnc_plot(tmp, CLIM, cellstr(num2str(comps)), (1:length(comps)), H(ii).H, 'Correlations (z)', sh(1), network_values, network_names);
            title(titleStr, 'parent', sh(1), 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            axis square;
            set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            
            if (nargout == 1)
                outFile = [dfncInfo.prefix, '_clusters_', icatb_returnFileIndex(ii), '.jpg'];
                printFile(H(ii).H, fullfile(html_dir, outFile));
                results(ii).file = outFile;
                results(ii).text = ['Number of occurrences of centroid ', num2str(ii), ' is ', titleStr];
            end
            
        end
        
%         %% Cluster Occurrences
%         
%         H(end + 1).H = icatb_getGraphics('Cluster Occurrences', 'graphics', 'dfnc_summary4', figVisible);
%         colormap(trd_fnc_cm_jetwhite([-1,1]));
%         
%         numCols = ceil(sqrt(dfncInfo.postprocess.num_clusters));
%         numRows = ceil(dfncInfo.postprocess.num_clusters/numCols);
%         
%         for ii = 1:dfncInfo.postprocess.num_clusters
%             sh = subplot(numRows, numCols, ii);
%             timeline = 0:(Nwin-1); timeline = timeline + dfncInfo.wsize/2; timeline = timeline*TR;
%             
%             for bb = 1:100
%                 pickME = ceil(M*rand(M,1));
%                 octemp = 100*mean(aIND(:, pickME) == ii,2);
%                 plot(timeline, octemp, 'Color', [.8 .8 .8])
%                 hold on
%             end
%             
%             hold on
%             oc = 100*mean(aIND == ii,2);
%             plot(timeline, oc, 'b')
%             if (ii == 1)
%                 xlabel('time (s)', 'parent', sh); ylabel('Number Of Occurrences (100 bootstraps)', 'parent', sh);
%             end
%             box off; set(sh, 'TickDir', 'out')
%             title(sprintf('%d (%d%%)', sum(clusterInfo.IDXall == ii),  round(100*sum(clusterInfo.IDXall == ii)/length(clusterInfo.IDXall))), 'horizontalAlignment', 'center', 'parent', sh, ...
%                 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
%             set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
%             
%         end
%         
%         if (nargout == 1)
%             outFile = [dfncInfo.prefix, '_cluster_occurrences.jpg'];
%             printFile(H(end).H, fullfile(html_dir, outFile));
%             results(end + 1).file = outFile;
%             results(end).text = 'No. of cluster occurrences using 100 bootstrap iterations';
%         end
        
        %% State vector stats
        [matNumTransitions, matFractStateTime, matTransitions, matMeanDwellTime] = m_state_vec_stats('_display_sum_scores.mat', dfncInfo.userInput.numOfSub, aIND, dfncInfo.postprocess.num_clusters, dfncInfo);

        H(end+1).H = icatb_getGraphics('State Vector Stats', 'graphics', 'dfnc_summary4', figVisible);
        colormap(jet);
        
        icatb_dfnc_plot_statevector_stats(dfncInfo.postprocess.num_clusters, matFractStateTime, matTransitions, matMeanDwellTime, matNumTransitions, H(end).H);
        
        set(findobj(H(end).H, 'type', 'axes'), 'YColor', FONT_COLOR, 'XColor', FONT_COLOR);
        
        try
            set(findobj(H(end).H, 'type', 'colorbar'), 'color', FONT_COLOR);
        catch
        end
        
        if (nargout == 1)
            
            %else
            titleStr = 'State Vector Stats';
            outFile = [dfncInfo.prefix, '_state_vector_stats.jpg'];
            printFile(H(end).H, fullfile(html_dir, outFile));
            results(end + 1).file = outFile;
            results(end).text = 'Plot shows frequency of each cluster followed by mean dwell time in windows followed by mean of state transition matrix across subjects';
        end
        
        %% State vector
        if (~isfield(clusterInfo, 'states'))
            % subjects x sessions x windows
            states = reshape(clusterInfo.IDXall, dfncInfo.userInput.numOfSess, dfncInfo.userInput.numOfSub, Nwin);
            states = permute(states, [2, 1, 3]);
        else
            states = clusterInfo.states;
            states = reshape(states, dfncInfo.userInput.numOfSub, dfncInfo.userInput.numOfSess, Nwin);
        end
        
        
%         numStateRows = ceil(sqrt(dfncInfo.userInput.numOfSub));
%         numStateCols = ceil(dfncInfo.userInput.numOfSub/numStateRows);
%         for nSess = 1:dfncInfo.userInput.numOfSess
%             H(end + 1).H = icatb_getGraphics(['State vector for session ', num2str(nSess)], 'graphics', ['dfnc_summary', num2str(4 + nSess)], figVisible);
%             for nSub = 1:dfncInfo.userInput.numOfSub
%                 sh = subplot(numStateRows, numStateCols, nSub);
%                 plot(squeeze(states(nSub, nSess, :)), 'm', 'parent', sh);
%                 title(['Sub ', num2str(nSub)], 'parent', sh, 'fontname', UI_FONTNAME, 'fontsize', 7);
%                 axis(sh, [0, Nwin, 0, dfncInfo.postprocess.num_clusters+1]);
%                 if (nSub ~= 1)
%                     set(sh, 'XTick', []);
%                     set(sh, 'YTick', []);
%                 end
%             end
%             set(findobj(H(end).H, 'type', 'axes'), 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
%             if (nargout == 1)
%                 titleStr = ['States for session ', num2str(nSess)];
%                 outFile = [dfncInfo.prefix, '_state_vector_session', num2str(nSess), '.jpg'];
%                 printFile(H(end).H, fullfile(html_dir, outFile));
%                 results(end + 1).file = outFile;
%                 results(end).text = ['State vector for session ', num2str(nSess)];
%             end
%         end    
    
    
    end



    if (nargout == 1)
        for nH = 1:length(H)
            close(H(nH).H);
        end
        varargout{1} = results;
        return;
    end

    icatb_plotNextPreviousExitButtons(H);
end


function printFile(H, file_name)

    pos = get(0, 'defaultFigurePosition');
    set(H, 'color', 'w');
    set(H, 'position', pos);
    set(H, 'resize', 'on');
    set(H, 'visible', 'on');
    titleH = get(findobj(H, 'type', 'axes'), 'title');
    if (iscell(titleH))
        % titleH = cell2mat(titleH);
        for nT = 1:length(titleH)
            set(titleH{nT}, 'color', 'k');
        end
    else
        set(titleH, 'color', 'k');
    end
    axesH = findobj(H, 'type', 'axes');
    set(axesH, 'YColor', 'k', 'XColor', 'k');
    set(axesH, 'fontname', 'times');
    set(axesH, 'fontsize', 8);
    textH = findobj(H, 'type', 'text');
    set(textH, 'color', 'k');
    set(textH, 'fontname', 'times');
    set(textH, 'fontsize', 8);
    try
        C = findobj(H, 'type', 'colorbar');
        set(C, 'color', 'k');
        set(get(C, 'label'), 'color', 'k');
    catch
    end
    %set(H,'PaperUnits','inches','PaperPosition', [0 0 4 3]);
    saveas(H, strrep(file_name, '.jpg', '.fig'));
    saveas(H, file_name);
end

function tc = regress_cov(tc, X)
    %% Regress covariates from corrs
    %

    betas = pinv(X)*tc;

    % Remove variance associated with the covariates
    tc = tc - X*betas;
end

function [p_masked, p]  = get_sig_pvalues(p, thresh, criteria)
    % apply fdr correction

    p_masked = thresh;
    if (strcmpi(criteria, 'mafdr'))
        p = mafdr(p);
    elseif(strcmpi(criteria, 'fdr') || strcmpi(criteria, 'bhfdr'))
        p_masked = icatb_fdr(p, thresh);
    end
end

function [matNumTransitions, matFractStateTime, matTransitions, matMeanDwellTime] = m_state_vec_stats(s_suffix, M, aIND, n_clusters, dfncInfo)
%STATS Summary of this function goes here
%   Detailed explanation goes here
% M is number of subjects
% aIND is the array holding the states for each window for all subjects
% dfncInfo are the variables from the dfnc parameter file


    matFractStateTime = zeros(M, n_clusters);
    matTransitions = zeros(M, n_clusters, n_clusters);
    matMeanDwellTime = zeros(M, n_clusters);
    matNumTransitions = zeros(M, 1);
    for ii = 1:M
        [FRii, TMii, MDTii, NTii] = icatb_dfnc_statevector_stats(aIND(:,ii), n_clusters);
        matFractStateTime(ii,:) = FRii;
        matTransitions(ii,:,:) = TMii;
        matMeanDwellTime(ii,:) = MDTii;
        matNumTransitions(ii) = NTii;
    end
    
    % save the summary statistics
    README_icatb={'Score Summary Information'; ...
            'Type following commands to get statistacal reports'; ...
            'mean(matFractStateTime(:,:)) % returns the average fractal state time for the group'; ...
            'squeeze(mean(matTransitions(:,:,:))) % returns the mean of the transition matrix' ; ...
            'mean(matMeanDwellTime) % returns the mean of the subject mean dwell times' ; ...
            'matNumTransitions % returns the number of transitions per subject'};
    save(fullfile(dfncInfo.outputDir, [dfncInfo.prefix s_suffix]), 'README_icatb', 'matNumTransitions', 'matFractStateTime', 'matTransitions', 'matMeanDwellTime');
    disp(['dFNC summary statistics saved to ' fullfile(dfncInfo.outputDir, [dfncInfo.prefix s_suffix])]);
end