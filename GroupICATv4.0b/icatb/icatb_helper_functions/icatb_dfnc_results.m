function varargout = icatb_dfnc_results(dfncInfo, display_criteria, statsInfo)
%% DFNC results
%

icatb_defaults;
global FONT_COLOR;
global UI_FONTNAME;
global UI_FS;

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

if (strcmpi(display_criteria, 'fnc oscillations'))
    
    %% FNC Oscillations
    H(1).H = icatb_getGraphics('FNC Oscillations (1)', 'graphics', 'dfnc_summary', figVisible);
    
    colormap(jet);
    
    %sh = subplot(1, 1, 1);
    sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
    
    CLIM = [min(FNCamp(:)), max(FNCamp(:))];
    FNCamp = icatb_vec2mat(FNCamp, 1);
    icatb_plot_FNC(FNCamp, CLIM, cellstr(num2str(comps)), (1:length(comps)), H(1).H, 'Average Std Of FNC Spectra', ...
        sh, network_values, network_names);
    title('Average Std Of FNC Spectra', 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
    set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
    
    if (nargout == 1)
        outFile = [dfncInfo.prefix, '_avg_std_oscillations.jpg'];
        printFile(H(1).H, fullfile(html_dir, outFile));
        results(1).file = outFile;
        results(1).text = 'Spectra is computed on dFNC correlations. Standard deviation of spectra is then averaged across subjects.';
    end
    
    
    H(2).H = icatb_getGraphics('FNC Oscillations (2)', 'graphics', 'dfnc_summary', figVisible);
    
    colormap(jet);
    
    sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
    
    CLIM = [min(FNCcm(:)), max(FNCcm(:))];
    FNCcm = icatb_vec2mat(FNCcm, 1);
    icatb_plot_FNC(FNCcm, CLIM, cellstr(num2str(comps)), (1:length(comps)), H(2).H, 'Average Center-Of-Mass Of FNC Spectra', sh, ...
        network_values, network_names);
    title('Average Center-Of-Mass Of FNC Spectra', 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
    set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
    
    if (nargout == 1)
        outFile = [dfncInfo.prefix, '_avg_center_of_mass_oscillations.jpg'];
        printFile(H(end).H, fullfile(html_dir, outFile));
        results(2).file = outFile;
        results(2).text = 'Spectral center of mass is computed on dFNC correlations and averaged across subjects.';
    end
    
elseif (strcmpi(display_criteria, 'group comparisons'))
    %% Group comparisons
    
    figVisible = 'on';
    
    resultFiles = {[dfncInfo.prefix, '_one_sample_ttest_results.mat'], [dfncInfo.prefix, '_paired_ttest_results.mat'], [dfncInfo.prefix, '_two_sample_ttest_results.mat']};
    
    countR = 0;
    for nR = 1:length(resultFiles)
        
        fname = fullfile(cluster_stats_directory, resultFiles{nR});
        if (exist(fname, 'file'))
            
            load(fname, 'mean_u', 't_u', 'p_u', 'groupNames', 'N', 'groupVals');
            
            %% Mean
            CLIM = abs([mean_u{:}]);
            CLIM = max(abs(CLIM(:)));
            CLIM = icatb_z_to_r(CLIM);
            CLIM = [-CLIM, CLIM];
            
            %% Loop over number of groups
            for nRows = 1:size(mean_u, 1)
                files = cell(1, size(mean_u, 2));
                %% Loop over number of cluster states
                for nC = 1:size(mean_u, 2)
                    tmpM = mean_u{nRows, nC};
                    if (~isempty(tmpM))
                        tmpM = reshape(icatb_z_to_r(tmpM(:)), size(tmpM));
                        H = icatb_getGraphics('Cluster stats', 'graphics', 'dfnc_summary', figVisible);
                        colormap(jet);
                        set(H, 'resize', 'on');
                        sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
                        tmpM = icatb_vec2mat(tmpM, 1);
                        titleText = ['State ', num2str(nC), ' ', groupNames{nRows}, '=', num2str(N(nRows, nC))];
                        [FH,AH,CH] = icatb_plot_FNC(tmpM, CLIM, cellstr(num2str(comps)), (1:length(comps)), H, titleText, sh, ...
                            network_values, network_names);
                        ylabel(CH, 'Correlations');
                        title(titleText, 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
                        set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
                        % write out mean maps
                        if strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_one_sample_ttest_results.mat'])
                            outFile = [dfncInfo.prefix, '_mean_cluster_ttest_', icatb_returnFileIndex(nC), '.jpg'];
                            tag = 'Mean correlations';
                        elseif strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_two_sample_ttest_results.mat'])
                            outFile = [dfncInfo.prefix, '_mean_cluster_ttest2_group', num2str(nRows), '_', icatb_returnFileIndex(nC), '.jpg'];
                            tag = ['Mean correlations of group ',  groupNames{nRows}];
                        else
                            outFile = [dfncInfo.prefix, '_mean_cluster_paired_ttest_cond', num2str(nRows), '_', icatb_returnFileIndex(nC), '.jpg'];
                            tag = ['Mean correlations of condition ',  groupNames{nRows}];
                        end
                        
                        printFile(H, fullfile(html_dir, outFile));
                        files{nC} = outFile;
                        close(H);
                        %end
                    end
                end
                %% End of loop over number of cluster states
                countR = countR + 1;
                results(countR).file = files;
                results(countR).text = ['Mean is computed for each cluster state of group ', groupNames{nRows}, '. Number of subjects with finite correlations is also shown.'];
                results(countR).tag = tag;
                clear tag;
            end
            %% End of loop over number of groups
            
            %% T-value
            files = cell(1, length(t_u));
            for nC = 1:length(t_u)
                %% use -1og10(p)*sign(t)
                ps = p_u{nC};
                ts = t_u{nC};
                if (~isempty(ts))
                    if (strcmpi(statsInfo.threshdesc, 'fdr'))
                        p_masked = icatb_fdr(ps, statsInfo.p_threshold);
                    else
                        p_masked = statsInfo.p_threshold;
                    end
                    good_inds = ps < p_masked;
                    
                    if (~isempty(find(good_inds == 1)))
                        
                        tmpM = zeros(size(ps));
                        tmpM(good_inds) = -log10(ps(good_inds) + eps).*sign(ts(good_inds));
                        tmpM = icatb_vec2mat(tmpM, 1);
                        tmpM(isnan(tmpM)) = 0;
                        
                        H = icatb_getGraphics('Cluster stats', 'graphics', 'dfnc_summary', figVisible);
                        colormap(jet);
                        set(H, 'resize', 'on');
                        sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
                        if (strcmpi(statsInfo.threshdesc, 'fdr'))
                            p_masked  = icatb_fdr(ps(find(isnan(ps) == 0)), statsInfo.p_threshold);
                            disp(p_masked)
                            fdrlim = -log10(p_masked);
                        else
                            p_masked = statsInfo.p_threshold;
                        end
                        
                        CLIM = [-max(abs(tmpM(:))) max(abs(tmpM(:)))];
                        titleText = ['Results of State ', num2str(nC)];
                        [FH,AH,CH] = icatb_plot_FNC(tmpM, CLIM, cellstr(num2str(comps)), (1:length(comps)), H, titleText, sh, ...
                            network_values, network_names);
                        ytick = get(CH, 'YTick');
                        if (strcmpi(statsInfo.threshdesc, 'fdr'))
                            set(CH, 'YTick', sort([ytick, -fdrlim fdrlim]));
                        end
                        ylabel(CH, '-sign(t) log_1_0(p-value)', 'Interpreter', 'tex');
                        title(titleText, 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
                        set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
                        % write out t-values
                        if strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_one_sample_ttest_results.mat'])
                            outFile = [dfncInfo.prefix, '_logp_value_', num2str(statsInfo.p_threshold), '_', statsInfo.threshdesc, '_cluster_ttest_', icatb_returnFileIndex(nC), '.jpg'];
                        elseif strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_two_sample_ttest_results.mat'])
                            outFile = [dfncInfo.prefix, '_logp_value_', num2str(statsInfo.p_threshold), '_', statsInfo.threshdesc, '_cluster_ttest2_', icatb_returnFileIndex(nC), '.jpg'];
                        else
                            outFile = [dfncInfo.prefix, '_logp_value_', num2str(statsInfo.p_threshold), '_', statsInfo.threshdesc, '_cluster_paired_ttest_', icatb_returnFileIndex(nC), '.jpg'];
                        end
                        printFile(H, fullfile(html_dir, outFile));
                        files{nC} = outFile;
                        close(H);
                    end
                end
            end
            
            if strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_one_sample_ttest_results.mat'])
                titleText = ['One sample t-test results of group ', groupNames{1}, ':'];
                tag = 'One sample t-test results';
            elseif strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_two_sample_ttest_results.mat'])
                titleText = ['Two sample t-test results of group ', groupNames{1}, ' vs ', ' group ', groupNames{2}, ':'];
                tag = 'Two sample t-test results';
            else
                titleText = ['Paired t-test results of condition ', groupNames{1}, ' vs ', ' condition ', groupNames{2}, ':'];
                tag = 'Paired t-test results';
            end
            
            countR = countR + 1;
            results(countR).file = files;
            results(countR).text = titleText;
            results(countR).tag = tag;
            clear tag;
            
            
            %% Plot mean dwell time for each group
            cluster_stats_file = fullfile(cluster_stats_directory,  [dfncInfo.prefix, '_cluster_stats.mat']);
            load(cluster_stats_file, 'state_vector_stats');
            mean_dwell_time = reshape (state_vector_stats.mean_dwell_time, dfncInfo.userInput.numOfSess, dfncInfo.userInput.numOfSub, size(state_vector_stats.mean_dwell_time, 2));
            H = icatb_getGraphics('Cluster stats', 'graphics', 'dfnc_summary', figVisible);
            colormap(jet);
            set(H, 'resize', 'on');
            sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
            colors = {'r', 'k'};
            area_colors = {[1, 0.6, 0.78], [0.8, 0.8, 0.8]};
            numClusterStates = size(state_vector_stats.mean_dwell_time, 2);
            if (size(mean_dwell_time, 1) > 1)
                % average across sessions
                mean_dwell_time = mean(mean_dwell_time);
            end
            mean_dwell_time = squeeze(mean_dwell_time);
            legendStr = cell(1, 2*length(groupVals));
            for nGroups = 1:length(groupVals)
                icatb_plot_with_ste_area(sh, 1:numClusterStates, mean_dwell_time(groupVals{nGroups}, :), [], colors{nGroups}, area_colors{nGroups});
                legendStr{2*nGroups - 1} = [groupNames{nGroups}, '+/-SEM'];
                legendStr{2*nGroups} = groupNames{nGroups};
                hold on;
            end
            legend(legendStr{:});
            hold off;
            title('Mean Dwell Time Vs Cluster States', 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            titleText = 'Mean dwell time in windows';
            ylabel(titleText);
            xlabel('Cluster states');
            if strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_one_sample_ttest_results.mat'])
                outFile = [dfncInfo.prefix, '_mdt_ttest.jpg'];
                tag = 'Mean dwell time';
            elseif strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_two_sample_ttest_results.mat'])
                outFile = [dfncInfo.prefix, '_mdt_ttest2.jpg'];
                tag = 'Mean dwell time of each group';
            else
                outFile = [dfncInfo.prefix, '_mdt_paired_ttest.jpg'];
                tag = 'Mean dwell time of each condition';
            end
            printFile(H, fullfile(html_dir, outFile));
            clear files;
            files{1} = outFile;
            countR = countR + 1;
            results(countR).file = files;
            results(countR).text = titleText;
            results(countR).tag = tag;
            clear tag;
            close(H);
            
            
            clear files;
            if strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_two_sample_ttest_results.mat'])
                allCombs = nchoosek(1:length(groupVals), 2);
                for nComb = 1:size(allCombs, 1)
                    
                    g1Name = groupNames{allCombs(nComb, 1)};
                    g2Name = groupNames{allCombs(nComb, 2)};
                    mG1 = mean_dwell_time(groupVals{allCombs(nComb, 1)}, :);
                    mG2 = mean_dwell_time(groupVals{allCombs(nComb, 2)}, :);
                    
                    pvals_mdt = zeros(1, size(mean_dwell_time, 2));
                    tvals_mdt = zeros(1, size(mean_dwell_time, 2));
                    
                    for nState = 1:length(pvals_mdt)
                        [pvals_mdt(nState), tvals_mdt(nState)] = icatb_ttest2(mG1(:, nState), mG2(:, nState), 0);
                    end
                    
                    titleStr = ['Two sample t-test of mean dwell time between groups ', g1Name, ' & ', g2Name];
                    
                    g1Name(isspace(g1Name)) = '';
                    g2Name(isspace(g2Name)) = '';
                    
                    ttest2FileName = [dfncInfo.prefix, '_ttest2_', g1Name, '_', g2Name, '.txt'];
                    numPara = 1;
                    varStruct(numPara).tag = 'State#';
                    varStruct(numPara).value = (1:length(pvals_mdt))';
                    
                    numPara = numPara + 1;
                    varStruct(numPara).tag = 'p-value';
                    varStruct(numPara).value = num2str(pvals_mdt(:), '%0.4f');
                    
                    numPara = numPara + 1;
                    varStruct(numPara).tag = 'T-value';
                    varStruct(numPara).value = num2str(tvals_mdt(:), '%0.4f');
                    
                    icatb_printToFile(fullfile(html_dir,  ttest2FileName), varStruct, '', 'column_wise');
                    clear varStruct;
                    
                    tag = ['Two sample t-test on MDT between ', g1Name, ' vs ', g2Name];
                    countR = countR + 1;
                    results(countR).file = {ttest2FileName};
                    results(countR).text = titleStr;
                    results(countR).tag = tag;
                    clear tag;
                end
            end
            
            
        end
        
    end
    
    if (~exist('results', 'var'))
        error('No results to display');
    end
    
    
    if (exist('meta_states_info', 'var'))
        metaFields = {'Num States', 'Change Between States', 'State Span', 'Total distance'};
        fNs_meta = fieldnames(meta_states_info);
        metaMethods = fNs_meta;
        metaMethods = strrep(metaMethods, 'pca', 'PCA');
        metaMethods = strrep(metaMethods, 'tica', 'Temporal ICA');
        metaMethods = strrep(metaMethods, 'sica', 'Spatial ICA');
        metaMethods = strrep(metaMethods, 'kmeans', 'K-means');
        dat = cell(1, length(fNs_meta));
        for nMeta = 1:length(fNs_meta)
            tmp =  [meta_states_info.(fNs_meta{nMeta}).summary.num_states, meta_states_info.(fNs_meta{nMeta}).summary.change_states, meta_states_info.(fNs_meta{nMeta}).summary.state_span, ...
                meta_states_info.(fNs_meta{nMeta}).summary.total_dist];
            dat{nMeta} = tmp;
            clear tmp;
        end
        %metaMethods = {'PCA', 'Temporal ICA', 'Spatial ICA', 'K-means'};
        % pca_dat = [meta_states_info.pca.summary.num_states, meta_states_info.pca.summary.change_states, meta_states_info.pca.summary.state_span, ...
        %     meta_states_info.pca.summary.total_dist];
        % tica_dat = [meta_states_info.tica.summary.num_states, meta_states_info.tica.summary.change_states, meta_states_info.tica.summary.state_span, ...
        %     meta_states_info.tica.summary.total_dist];
        % sica_dat = [meta_states_info.sica.summary.num_states, meta_states_info.sica.summary.change_states, meta_states_info.sica.summary.state_span, ...
        %     meta_states_info.sica.summary.total_dist];
        %kmeans_dat = [meta_states_info.kmeans.summary.num_states, meta_states_info.kmeans.summary.change_states, meta_states_info.kmeans.summary.state_span, ...
        %    meta_states_info.kmeans.summary.total_dist];
        %dat = {pca_dat, tica_dat, sica_dat, kmeans_dat};
        %clear pca_dat tica_dat sica_dat kmeans_dat;
        for nR = 1:length(resultFiles)
            fname = fullfile(cluster_stats_directory, resultFiles{nR});
            if (exist(fname, 'file'))
                countR = length(results) + 1;
                load(fname, 'groupNames', 'groupVals');
                if strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_one_sample_ttest_results.mat'])
                    
                    tag = 'Meta State: One sample t-test';
                    titleStr = ['Meta State: One sample t-test on ', groupNames{1}];
                    outFile = [dfncInfo.prefix, '_meta_state_ttest_', groupNames{1}, '.txt'];
                    
                    T = [];
                    P = [];
                    M = [];
                    clear t_u p_u stats_u;
                    for nM = 1:length(metaMethods)
                        tmpp1 = dat{nM}(groupVals{1}, :);
                        % [ t_u, p_u, stats_u] = mT(tmpp1, ones(length(groupVals{1}), 1), [], 0, {'verbose'});
                        t_u = zeros(1, size(tmpp1, 2));
                        p_u = t_u;
                        for nTcols = 1:size(tmpp1, 2)
                            [p_u(nTcols), t_u(nTcols)] = icatb_ttest(tmpp1(:, nTcols));
                            %[ t_u, p_u, stats_u] = mT(tmpp1, ones(length(groupVals{1}), 1), [], 0, {'verbose'});
                        end
                        T = [T; t_u];
                        P = [P; p_u];
                        M = [M; mean(tmpp1)];
                    end
                    
                    clear varStruct;
                    
                    numPara = 1;
                    varStruct(numPara).tag = 'Method';
                    varStruct(numPara).value = char(metaMethods);
                    
                    for nP = 1:size(T, 2)
                        numPara = numPara + 1;
                        varStruct(numPara).tag = metaFields{nP};
                        varStruct(numPara).value = strcat('Tval = ', num2str(T(:, nP),'%0.4f'), ', Pval = ', num2str(P(:, nP),'%0.4f'), ...
                            ', Mean = ', num2str(M(:, nP),'%0.4f'));
                    end
                    
                    icatb_printToFile(fullfile(html_dir,  outFile), varStruct, '', 'column_wise');
                    clear varStruct;
                elseif strcmpi(resultFiles{nR}, [dfncInfo.prefix, '_two_sample_ttest_results.mat'])
                    
                    outFile = [dfncInfo.prefix, '_meta_state_ttest2_', groupNames{1}, '_', groupNames{2}, '.txt'];
                    tag = 'Meta State: Two sample t-test';
                    titleStr = ['Meta State: Two sample t-test between groups ', groupNames{1},  ' and ', groupNames{2}];
                    
                    N1 = length(groupVals{1});
                    N2 = length(groupVals{2});
                    X = ones(N1 + N2, 1);
                    X(N1 + 1:end) = 0;
                    
                    T = [];
                    P = [];
                    M1 = [];
                    M2 = [];
                    clear t_u p_u stats_u;
                    
                    for nM = 1:length(metaMethods)
                        tmp1 = dat{nM}(groupVals{1}, :);
                        tmp2 = dat{nM}(groupVals{2}, :);
                        %[ t_u, p_u, stats_u] = mT(dat{nM}(groupVals{1}, :), X, [], 0, {'verbose'});
                        [t_u, p_u, stats_u] = mT([tmp1;tmp2],X, [], 1, {'verbose'});
                        T = [T; t_u];
                        P = [P; p_u];
                        M1 = [M1; mean(tmp1)];
                        M2 = [M2; mean(tmp2)];
                    end
                    
                    clear varStruct;
                    
                    numPara = 1;
                    varStruct(numPara).tag = 'Method';
                    varStruct(numPara).value = char(metaMethods);
                    
                    for nP = 1:size(T, 2)
                        numPara = numPara + 1;
                        varStruct(numPara).tag = metaFields{nP};
                        varStruct(numPara).value = strcat('Tval = ', num2str(T(:, nP),'%0.4f'), ', Pval = ', num2str(P(:, nP),'%0.4f'), ...
                            ', Mean of Group1 = ', num2str(M1(:, nP),'%0.4f'), ', Mean of Group2 = ', num2str(M2(:, nP),'%0.4f'));
                    end
                    
                    icatb_printToFile(fullfile(html_dir,  outFile), varStruct, '', 'column_wise');
                    clear varStruct;
                else
                    
                    outFile = [dfncInfo.prefix, '_meta_state_ttest_', groupNames{1}, '_', groupNames{2}, '.txt'];
                    titleStr = ['Meta State: Paired t-test between conditions ', groupNames{1},  ' and ', groupNames{2}];
                    tag = 'Meta State: Paired t-test';
                    
                    N1 = length(groupVals{1});
                    N2 = length(groupVals{2});
                    X = ones(N1, 1);
                    
                    T = [];
                    P = [];
                    M1 = [];
                    M2 = [];
                    clear t_u p_u stats_u;
                    
                    for nM = 1:length(metaMethods)
                        tmp1 = dat{nM}(groupVals{1}, :);
                        tmp2 = dat{nM}(groupVals{2}, :);
                        [t_u, p_u, stats_u] = mT(tmp1 - tmp2, X, [], 0, {'verbose'});
                        T = [T; t_u];
                        P = [P; p_u];
                        M1 = [M1; mean(tmp1)];
                        M2 = [M2; mean(tmp2)];
                    end
                    
                    clear varStruct;
                    
                    numPara = 1;
                    varStruct(numPara).tag = 'Method';
                    varStruct(numPara).value = char(metaMethods);
                    
                    for nP = 1:size(T, 2)
                        numPara = numPara + 1;
                        varStruct(numPara).tag = metaFields{nP};
                        varStruct(numPara).value = strcat('Tval = ', num2str(T(:, nP),'%0.4f'), ', Pval = ', num2str(P(:, nP),'%0.4f'), ...
                            ', Mean of Cond1 = ', num2str(M1(:, nP),'%0.4f'), ', Mean of Cond2 = ', num2str(M2(:, nP),'%0.4f'));
                    end
                    
                    icatb_printToFile(fullfile(html_dir,  outFile), varStruct, '', 'column_wise');
                    clear varStruct;
                    
                    
                end
                
                results(countR).file = {outFile};
                results(countR).text = titleStr;
                results(countR).tag = tag;
                
            end
        end
    end
    
    varargout{1} = results;
    return;
    
    
    %H(1).H = icatb_getGraphics('FNC Oscillations (1)', 'graphics', 'dfnc_summary', figVisible);
    
    
elseif strcmpi(display_criteria, 'meta state analysis')
    
    if (~exist('dfncInfo', 'var'))
        error('Meta states info doesn''t exist. You need to re-run the post-processing step');
    end
    
    
    %% Meta states (pca, temporal ica, spatial ica)
    fNs_meta = fieldnames(meta_states_info);
    metaMethods = fNs_meta;
    metaMethods = strrep(metaMethods, 'pca', 'PCA');
    metaMethods = strrep(metaMethods, 'tica', 'Temporal ICA');
    metaMethods = strrep(metaMethods, 'sica', 'Spatial ICA');
    metaMethods = strrep(metaMethods, 'kmeans', 'K-means');
    cp_data = cell(1, length(fNs_meta));
    for nMeta = 1:length(fNs_meta)
        cp_data{nMeta} = meta_states_info.(fNs_meta{nMeta}).S;
    end
    
    numComps = size(cp_data{1}, 2);
    titles = metaMethods;
    %titles = {'PCA', 'Temporal ICA', 'Spatial ICA', 'K-means'};
    %cp_data = {meta_states_info.pca.S, meta_states_info.tica.S, meta_states_info.sica.S, meta_states_info.kmeans.S};
    %numComps = size(meta_states_info.sica.S, 2);
    %numRows = ceil(sqrt(numComps));
    %numCols = ceil(numComps/numRows);
    figCount = 0;
    for nFig = 1:length(titles)
        %countS = 0;
        %         H(nFig).H = icatb_getGraphics(['Correlation patterns (', titles{nFig} ')'], 'graphics', 'dfnc_summary', figVisible);
        %         colormap(jet);
        tmp = cp_data{nFig};
        CLIM = max(abs(tmp(:)));
        CLIM = [-CLIM, CLIM];
        
        for countS = 1:numComps
            figCount = figCount + 1;
            H(figCount).H = icatb_getGraphics(['Correlation patterns (', titles{nFig}, num2str(countS), ')'], 'graphics', 'dfnc_summary', figVisible);
            colormap(jet);
            tmp2 = icatb_vec2mat(tmp(:, countS), 1);
            sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
            if (strcmpi(titles{nFig}, 'k-means'))
                cbarTitle = 'Correlations (z)';
            else
                cbarTitle = 'Arbitrary';
            end
            [F, A, C, I] = icatb_plot_FNC(tmp2, CLIM, cellstr(num2str(comps)), (1:length(comps)), H(figCount).H, cbarTitle, ...
                sh, network_values, network_names);
            title([titles{nFig}, ' CP ', num2str(countS)], 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
            if (nargout == 1)
                outFile = [dfncInfo.prefix, '_meta_states_', num2str(nFig), '_state_', num2str(countS), '.jpg'];
                printFile(H(figCount).H, fullfile(html_dir, outFile));
                results(figCount).file = outFile;
                if (countS == 1)
                    results(figCount).text = [titles{nFig}, ' is computed on the windowed dFNC correlations (windows*subjects by component pairs).'];
                else
                    results(figCount).text = '';
                end
            end
            
        end
        
        %         for nRow = 1:numRows
        %             for nCol = 1:numCols
        %                 countS = countS + 1;
        %                 if (countS <= numComps)
        %                     sh = subplot(numRows, numCols, countS);
        %                     tmp2 = icatb_vec2mat(tmp(:, countS), 1);
        %                     [F, A, C, I] = icatb_plot_FNC(tmp2, CLIM, cellstr(num2str(comps)), (1:length(comps)), H(nFig).H, 'Arbitrary', ...
        %                         sh, network_values, network_names);
        %                     if (countS ~= numComps)
        %                         delete(C);
        %                     else
        %                         cpos = get(C, 'position');
        %                         cpos(1) = cpos(1) + 0.05;
        %                         set(C, 'position', cpos);
        %                     end
        %                     clear F A C I;
        %                     title([titles{nFig}, ' CP ', num2str(countS)], 'parent', sh, 'horizontalAlignment', 'center', 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
        %                     set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
        %                 end
        %             end
        %         end
        
        %         if (nargout == 1)
        %             outFile = [dfncInfo.prefix, '_meta_states_', num2str(nFig), '.jpg'];
        %             printFile(H(nFig).H, fullfile(html_dir, outFile));
        %             results(nFig).file = outFile;
        %             results(nFig).text = [titles{nFig}, ' is computed on the windowed dFNC correlations (windows*subjects by component pairs).'];
        %         end
        
    end
    
else
    
    %% Cluster centroids
    M = length (dfncInfo.outputFiles);
    Nwin = length(clusterInfo.IDXall) / M;
    aIND = reshape(clusterInfo.IDXall, M, Nwin);
    aIND = aIND';
    for ii = 1:dfncInfo.postprocess.num_clusters
        fig_title = ['Cluster Centroid (', num2str(ii), ')'];
        H(ii).H = icatb_getGraphics(fig_title, 'graphics', 'dfnc_summary3', figVisible);
        colormap(jet);
        sh = axes('units', 'normalized', 'position', [0.2, 0.2, 0.6, 0.6]);
        CLIM = max(abs(clusterInfo.Call(:)));
        CLIM = [-CLIM, CLIM];
        tmp = icatb_vec2mat(clusterInfo.Call(ii, :), 1);
        titleStr = sprintf('%d (%d%%)', sum(clusterInfo.IDXall == ii),  round(100*sum(clusterInfo.IDXall == ii)/length(clusterInfo.IDXall)));
        icatb_plot_FNC(tmp, CLIM, cellstr(num2str(comps)), (1:length(comps)), H(ii).H, 'Correlations (z)', sh(1), network_values, network_names);
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
    
    %% Cluster Occurrences
    
    H(end + 1).H = icatb_getGraphics('Cluster Occurrences', 'graphics', 'dfnc_summary4', figVisible);
    colormap(jet);
    
    numCols = ceil(sqrt(dfncInfo.postprocess.num_clusters));
    numRows = ceil(dfncInfo.postprocess.num_clusters/numCols);
    
    for ii = 1:dfncInfo.postprocess.num_clusters
        sh = subplot(numRows, numCols, ii);
        timeline = 0:(Nwin-1); timeline = timeline + dfncInfo.wsize/2; timeline = timeline*TR;
        
        for bb = 1:100
            pickME = ceil(M*rand(M,1));
            octemp = 100*mean(aIND(:, pickME) == ii,2);
            plot(timeline, octemp, 'Color', [.8 .8 .8])
            hold on
        end
        
        hold on
        oc = 100*mean(aIND == ii,2);
        plot(timeline, oc, 'b')
        if (ii == 1)
            xlabel('time (s)', 'parent', sh); ylabel('Number Of Occurrences (100 bootstraps)', 'parent', sh);
        end
        box off; set(sh, 'TickDir', 'out')
        title(sprintf('%d (%d%%)', sum(clusterInfo.IDXall == ii),  round(100*sum(clusterInfo.IDXall == ii)/length(clusterInfo.IDXall))), 'horizontalAlignment', 'center', 'parent', sh, ...
            'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
        set(sh, 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
        
    end
    
    if (nargout == 1)
        outFile = [dfncInfo.prefix, '_cluster_occurrences.jpg'];
        printFile(H(end).H, fullfile(html_dir, outFile));
        results(end + 1).file = outFile;
        results(end).text = 'No. of cluster occurrences using 100 bootstrap iterations';
    end
    
    %% State vector stats
    
    aFR = zeros(M, dfncInfo.postprocess.num_clusters);
    aTM = zeros(M, dfncInfo.postprocess.num_clusters, dfncInfo.postprocess.num_clusters);
    aMDT = zeros(M, dfncInfo.postprocess.num_clusters);
    aNT = zeros(M, 1);
    for ii = 1:M
        [FRii, TMii, MDTii, NTii] = icatb_dfnc_statevector_stats(aIND(:,ii), dfncInfo.postprocess.num_clusters);
        aFR(ii,:) = FRii;
        aTM(ii,:,:) = TMii;
        aMDT(ii,:) = MDTii;
        aNT(ii) = NTii;
    end
    
    H(end+1).H = icatb_getGraphics('State Vector Stats', 'graphics', 'dfnc_summary4', figVisible);
    colormap(jet);
    
    icatb_dfnc_plot_statevector_stats(dfncInfo.postprocess.num_clusters, aFR, aTM, aMDT, aNT, H(end).H);
    
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
    
    
    numStateRows = ceil(sqrt(dfncInfo.userInput.numOfSub));
    numStateCols = ceil(dfncInfo.userInput.numOfSub/numStateRows);
    for nSess = 1:dfncInfo.userInput.numOfSess
        H(end + 1).H = icatb_getGraphics(['State vector for session ', num2str(nSess)], 'graphics', ['dfnc_summary', num2str(4 + nSess)], figVisible);
        for nSub = 1:dfncInfo.userInput.numOfSub
            sh = subplot(numStateRows, numStateCols, nSub);
            plot(squeeze(states(nSub, nSess, :)), 'm', 'parent', sh);
            title(['Sub ', num2str(nSub)], 'parent', sh, 'fontname', UI_FONTNAME, 'fontsize', 7);
            axis(sh, [0, Nwin, 0, dfncInfo.postprocess.num_clusters+1]);
            if (nSub ~= 1)
                set(sh, 'XTick', []);
                set(sh, 'YTick', []);
            end
        end
        set(findobj(H(end).H, 'type', 'axes'), 'YColor', FONT_COLOR, 'XColor', FONT_COLOR, 'fontname', UI_FONTNAME, 'fontsize', UI_FS-2);
        if (nargout == 1)
            titleStr = ['States for session ', num2str(nSess)];
            outFile = [dfncInfo.prefix, '_state_vector_session', num2str(nSess), '.jpg'];
            printFile(H(end).H, fullfile(html_dir, outFile));
            results(end + 1).file = outFile;
            results(end).text = ['State vector for session ', num2str(nSess)];
        end
    end
    
    
end



if (nargout == 1)
    for nH = 1:length(H)
        close(H(nH).H);
    end
    varargout{1} = results;
    return;
end

icatb_plotNextPreviousExitButtons(H);


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
