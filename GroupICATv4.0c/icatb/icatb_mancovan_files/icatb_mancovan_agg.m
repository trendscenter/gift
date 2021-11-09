function icatb_mancovan_agg(mancovanInfo)
%% Mancova aggregate function. Saves feature and univariate stats output to be used in later in mancova reduce function
%

disp('Storing feature and stats info in *stats*info* file ... ');
outputDir = mancovanInfo.outputDir;

try
    freq_limits = mancovanInfo.display.freq_limits;
catch
    freq_limits = [0.1, 0.15];
end

load (mancovanInfo.outputFiles(1).filesInfo(1).result_files{1}, 'UNI');
tests = UNI.tests;

feature_info = cell(length(mancovanInfo.outputFiles), length(mancovanInfo.comps));
uni_results_info = cell(length(mancovanInfo.outputFiles), length(mancovanInfo.comps), length(tests));

if (~isempty(strmatch('spatial maps', lower(mancovanInfo.features), 'exact')))
    if (strcmpi(mancovanInfo.userInput.feature_params.final.sm_mask, 'default') || isempty(mancovanInfo.userInput.feature_params.final.sm_mask))
        warning('You need to use explicit mask in oder to use distributed mancova. Otherwise, skip spatial maps feature and use rest of the features');
        return;
    end
end

for nF = 1:length(mancovanInfo.outputFiles)
    
    disp(['Working on feature ', mancovanInfo.features{nF}, ' ...']);
    result_files = mancovanInfo.outputFiles(nF).filesInfo(1).result_files;
    
    if (strcmpi(mancovanInfo.features{nF}, 'spatial maps'))
        varToLoad = 'SM';
    elseif (strcmpi(mancovanInfo.features{nF}, 'timecourses spectra'))
        varToLoad = 'spectra_tc';
    elseif (strcmpi(mancovanInfo.features{nF}, 'fnc correlations'))
        varToLoad = 'fnc_corrs';
    else
        varToLoad = 'fnc_values';
    end
    
    
    for nR = 1:length(result_files)
        tmp = load(fullfile(outputDir, result_files{nR}), varToLoad);
        tmp = tmp.(varToLoad);
        sum_feature = sum(tmp);
        ssq_features = sum(tmp.^2);
        load(fullfile(outputDir, result_files{nR}), 'UNI');
        out.sum_feature = sum_feature;
        out.ssq_features = ssq_features;
        out.N = size(tmp, 1);
        if (nR == 1)
            if (strcmpi(varToLoad, 'spectra_tc'))
                load(fullfile(outputDir, result_files{nR}), 'freq');
                % if (nR == 1)
                out.freq = freq;
                % end
            elseif  (strcmpi(varToLoad, 'SM'))
                load(fullfile(outputDir, result_files{nR}), 'mask_ind');
                out.mask_ind = mask_ind;
            end
        end
        
        varInfo = whos('-file', fullfile(outputDir, result_files{nR}));
        if (~isempty(strmatch('low_inds', cellstr(char(varInfo.name)), 'exact')))
            load(fullfile(outputDir, result_files{nR}), 'low_inds')
            out.low_inds = low_inds;
            clear low_inds;
        end
        
        load(fullfile(outputDir, result_files{nR}), 'comp_number');
        out.comp_number = comp_number;
        
        
        if (strcmpi(varToLoad, 'spectra_tc'))
            %% Spectra tc
            fALFF = 0;
            dynamicrange = 0;
            for nS = 1:size(tmp, 1)
                [tmp_dyn_range, tmp_falff] = icatb_get_spec_stats(tmp(nS, :), freq, freq_limits);
                fALFF = fALFF + tmp_falff;
                dynamicrange = dynamicrange + tmp_dyn_range;
            end
            out.dynamic_range = dynamicrange;
            out.fALFF = fALFF;
        end
        
        feature_info{nF, nR} = out;
        
        clear out;
        
        % Univariate results info
        for nTest = 1:length(tests)
            uni_results_info{nF, nR, nTest} = UNI.stats{nTest};
        end
        
        
    end
    
end

mancova_agg_file = fullfile(outputDir, [mancovanInfo.prefix, '_stats_info.mat']);

mInfo = mancovanInfo;
clear mancovanInfo;
mancovanInfo.comps = mInfo.comps;
mancovanInfo.features = mInfo.features;
mancovanInfo.outputFiles = mInfo.outputFiles;
mancovanInfo.cov = mInfo.cov;
mancovanInfo.prefix = mInfo.prefix;
mancovanInfo.comp = mInfo.comp;
mancovanInfo.tests = tests;
try
    mancovanInfo.HInfo = mInfo.HInfo;
catch
end

try
    mancovanInfo.time = mInfo.time;
catch
end

try
    mancovanInfo.univInfo = mInfo.univInfo;
catch
end


try
    mancovanInfo.display = mInfo.display;
catch
end

p_threshold =  mInfo.userInput.p_threshold;
try
    p_threshold = mInfo.display.p_threshold;
catch
end

mancovanInfo.display.p_threshold = p_threshold;

save(mancova_agg_file, 'feature_info', 'uni_results_info', 'mancovanInfo', '-v7.3');

disp('Done');
fprintf('\n');