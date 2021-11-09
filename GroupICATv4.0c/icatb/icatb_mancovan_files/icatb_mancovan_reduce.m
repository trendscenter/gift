function mancovanInfo = icatb_mancovan_reduce(files, outputDir, displayInfo)
%% Mancova reduce. Compute global output from local stats output
%
% Inputs:
% files - Mancova param files in a cell array
% outputDir - Output directory to place the gathered stats output
%

files = cellstr(files);

if (~exist('outputDir', 'var') || isempty(outputDir))
    outputDir = fullfile(pwd, 'tmp');
end

if (~exist(outputDir, 'dir'))
    mkdir(outputDir);
end

load(files{1}, 'mancovanInfo');
nfeatures = length(mancovanInfo.features);
mancovanInfo.outputDir = fileparts(files{1});
%numResults = length(mancovanInfo.comps);
%load (fullfile(mancovanInfo.outputDir, mancovanInfo.outputFiles(1).filesInfo(1).result_files{1}), 'UNI');
tests = mancovanInfo.tests;
nTests = length(tests);
%clear UNI;
mancovanInfo.outputDir = outputDir;
mancovanInfo.userInput.outputDir = outputDir;

if (exist('displayInfo', 'var') && ~isempty(displayInfo))
    mancovanInfo.display = displayInfo;
end

save(fullfile(outputDir, [mancovanInfo.prefix, '.mat']), 'mancovanInfo');

cov_names = cellstr(char(mancovanInfo.cov.name));
cov_labels = [mancovanInfo.cov.labels];

%% Aggregate statistics
% Loop over all features
for nF = 1:nfeatures
    
    numResults = length(mancovanInfo.outputFiles(nF).filesInfo(1).result_files);
    
    disp(['Computing global stats on feature ', mancovanInfo.features{nF}, ' ...']);
    
    % Loop over all components
    for nR = 1:numResults
        
        if (numResults > 1)
            comp_number = mancovanInfo.comps(nR);
        else
            comp_number = mancovanInfo.comps;
        end
        
        t = cell(1, nTests);
        p = cell(1, nTests);
        stats = cell(1, nTests);
        
        % Loop over all tests
        for nT = 1:nTests
            
            stats_all = cell(length(files), 1);
            % Loop over all sites
            for nFile = 1:length(files)
                fname = fullfile(fileparts(files{nFile}), [mancovanInfo.prefix, '_stats_info.mat']);
                matObj = matfile(fname);
                
                stats_all(nFile) = matObj.uni_results_info(nF, nR, nT);
            end
            
            % Call reduce function to gather global stats output from local stats
            [t{nT}, p{nT}, stats{nT}] = icatb_mT_reduce(stats_all, stats_all{1}.Terms, stats_all{1}.Term, {});
            addTag = 1;
            chkNames = strmatch(lower(tests{nT}), lower(cov_names), 'exact');
            if (~isempty(chkNames))
                if (strcmpi(mancovanInfo.cov(chkNames(1)).type, 'continuous'))
                    addTag = 0;
                end
            end
            stats{nT} = get_contrast_label(stats{nT}, tests{nT}, cov_names, cov_labels, addTag);
            
        end
        % End loop over all tests
        
        
        R2 = cell(size(t));
        for nU = 1:length(t)
            R2{nU} = t{nU}.^2 ./ (t{nU}.^2 + stats{nU}.DFE);
        end
        
        UNI.t = t;
        UNI.p = p;
        UNI.tests = tests;
        UNI.stats = stats;
        UNI.R2 = R2;
        MULT = [];
        
        chk = fileparts(mancovanInfo.outputFiles(nF).filesInfo(1).result_files{nR});
        if (~isempty(chk))
            if ~exist(fullfile(outputDir, chk), 'dir')
                mkdir(outputDir, chk);
            end
        end
        
        out_file = fullfile(outputDir, mancovanInfo.outputFiles(nF).filesInfo(1).result_files{nR});
        save(out_file, 'MULT', 'UNI', 'comp_number', '-v7.3');
        
        fname = fullfile(fileparts(files{1}), [mancovanInfo.prefix, '_stats_info.mat']);
        
        matObj = matfile(fname);
        
        if (strcmpi(mancovanInfo.features{nF}, 'spatial maps'))
            mask_ind = matObj.feature_info(nF, 1);
            mask_ind = mask_ind{1}.mask_ind;
            save(out_file, 'mask_ind', '-append');
        elseif (strcmpi(mancovanInfo.features{nF}, 'timecourses spectra'))
            
            freq = matObj.feature_info(nF, 1);
            freq = freq{1}.freq;
            save(out_file, 'freq', '-append');
        else
            if (nR == 2)
                linds = matObj.feature_info(nF, 2);
                try
                    low_inds = linds{1}.low_inds;
                    save(out_file, 'low_inds', '-append');
                catch
                end
            end
        end
        
    end
    % End loop over all components
end
% End Loop over all features



%% Aggregate feature summary
% Loop over features
for nF = 1:nfeatures
    
    disp(['Aggregating features ', mancovanInfo.features{nF}, ' ...']);
    
    numResults = length(mancovanInfo.outputFiles(nF).filesInfo(1).result_files);
    
    % Loop over all components
    for nR = 1:numResults
        
        num_subjects = 0;
        sum_feature = 0;
        ssq_features = 0;
        dynamic_range = 0;
        fALFF = 0;
        
        % Loop over files
        for nFile = 1:length(files)
            
            fname = fullfile(fileparts(files{nFile}), [mancovanInfo.prefix, '_stats_info.mat']);
            matObj = matfile(fname);
            out = matObj.feature_info(nF, nR);
            out = out{1};
            
            sum_feature = sum_feature + out.sum_feature;
            ssq_features = ssq_features + out.ssq_features;
            num_subjects = num_subjects + out.N;
            
            if (isfield(out, 'mask_ind'))
                mask_ind = out.mask_ind;
            end
            if (isfield(out, 'low_inds'))
                low_inds = out.low_inds;
            end
            
            if (strcmpi(mancovanInfo.features{nF}, 'timecourses spectra'))
                dynamic_range = dynamic_range + out.dynamic_range;
                fALFF = fALFF + out.fALFF;
            end
            
        end
        % End loop over files
        
        % compute mean and sem
        feature_mean = sum_feature / num_subjects;
        agg.feature.mean = feature_mean;
        agg.feature.sem = icatb_compute_sem(feature_mean, ssq_features, num_subjects);
        
        if (strcmpi(mancovanInfo.features{nF}, 'spatial maps'))
            % Write t-map files
            tmapVol = mancovanInfo.HInfo(1);
            tmapVol.n(1) = 1;
            tmap = zeros(mancovanInfo.HInfo.dim(1:3));
            tmap_file = fullfile(outputDir, mancovanInfo.outputFiles(nF).filesInfo.tmap_files{nR});
            tmap(mask_ind) =  agg.feature.mean ./  agg.feature.sem;
            tmapVol.fname = tmap_file;
            icatb_write_vol(tmapVol, tmap);
        elseif (strcmpi(mancovanInfo.features{nF}, 'timecourses spectra'))
            % append dynamic range and faLFF
            dynamic_range = dynamic_range / num_subjects;
            fALFF = fALFF / num_subjects;
            agg.feature.dynamic_range = dynamic_range;
            agg.feature.fALFF = fALFF;
        elseif (strcmpi(mancovanInfo.features{nF}, 'fnc correlations'))
            % FNC (lag and no-lag)
            if (nR == 1)
                FNCM = icatb_vec2mat(icatb_z_to_r(squeeze(agg.feature.mean)), 1);
            else
                % domain average FNCM
                FNCM = squeeze(icatb_low2Mat(icatb_z_to_r(squeeze(agg.feature.mean)), length(out.comp_number), low_inds));
            end
            agg.feature.mean = FNCM;
            
        else
            FNCM = icatb_vec2mat(icatb_z_to_r(squeeze(agg.feature.mean)), 1);
            agg.feature.mean = FNCM;
            
        end
        % end for adding aggregate fields
        
        out_file = fullfile(outputDir, mancovanInfo.outputFiles(nF).filesInfo(1).result_files{nR});
        save(out_file, 'agg', '-append');
        clear agg;
        
    end
    % End loop over components
end
% End loop over features

disp('Done');
fprintf('\n');

function s = get_contrast_label(s, tname, allnames, alllabels, addTag)

if (~exist('addTag', 'var'))
    addTag = 0;
end

prefixV = '';
if (addTag)
    prefixV = [tname, '_'];
end

if length(s.Term) == 1 %main effect
    varIND = find(strcmp(tname, allnames));
    labels = alllabels{varIND};
    if (~iscell(labels))
        labels = {labels};
    end
    for jj = 1:size(s.Levels,1)
        if s.Levels(jj,2) == 0 && length(labels) == 1
            s.Contrast{jj} = [prefixV, labels{1}];
        else
            s.Contrast{jj} = [prefixV, '(' labels{s.Levels(jj,1)+1} ') - (' labels{s.Levels(jj,2)+1} ')'];
            
        end
    end
    
else %interaction
    term1_end = strfind(tname, '_X_')-1;
    term2_start = term1_end + 4;
    clear varIND
    varIND(1) = find(strcmp(tname(1:term1_end), allnames));
    varIND(2) = find(strcmp(tname(term2_start:end), allnames));
    labels_1 = alllabels{varIND(1)};
    labels_2 = alllabels{varIND(2)};
    
    if (~iscell(labels_1))
        labels_1 = {labels_1};
    end
    
    if (~iscell(labels_2))
        labels_2 = {labels_2};
    end
    
    
    for jj = 1:size(s.Levels,1)
        if length(labels_1)*length(labels_2) == 1 % both are continuous variables
            s.Contrast{jj} = [prefixV, '(' labels_1{1} ') X (' labels_2{1} ')'];
        elseif length(labels_1) == 1 || length(labels_2) == 1 %one continuous, one categorical
            
            if length(labels_1) > length(labels_2)
                temp = labels_1;
                labels_1 = labels_2;
                labels_2 = temp;
            end
            
            s.Contrast{jj} = [prefixV, '(' labels_1{1} ') X [(' labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
        elseif length(labels_1)*length(labels_2) == 4 %categorical, two by two
            s.Contrast{jj} = [prefixV, '[(' labels_1{s.Levels(jj,1)+1} ') - (' labels_1{s.Levels(jj,2)+1} ')] X [(' ...
                labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
        elseif length(labels_1)*length(labels_2) == 6 %categorical, two by three
            if length(labels_1) > length(labels_2)
                temp = labels_1;
                labels_1 = labels_2;
                labels_2 = temp;
            end
            s.Contrast{jj} = [prefixV, '[(' labels_1{2} ') - (' labels_1{1} ')] X [(' labels_2{s.Levels(jj,1)+1} ') - (' labels_2{s.Levels(jj,2)+1} ')]'];
            
        else
            s.Contrast{jj} = tname;
            %             combinations = [];
            %             for j = 1 : length(labels_1)
            %                 for k = 1 :length(labels_2)
            %                     combinations(end + 1, :) = [ j k ];
            %                 end
            %             end
        end
    end
    
end
