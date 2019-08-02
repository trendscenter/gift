function icatb_postprocess_timecourses(param_file)
%% Write FNC and spectra information in *_postprocess_results.mat.
% FNC correlations (transformed to fisher z-scores) are saved as fnc_corrs_all with
% dimensions subjects x sessions x components x components. Spectra is
% saved as spectra_tc_all of dimensions subjects x sessions x spectral length.
%

%% Load defaults
icatb_defaults;
global EXPERIMENTAL_TR;
global TIMECOURSE_POSTPROCESS;
global DETRENDNUMBER;
global PARAMETER_INFO_MAT_FILE;

filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
if (~exist('param_file', 'var'))
    param_file = icatb_selectEntry('typeEntity', 'file', 'title', 'Select Parameter File', 'filter', filterP);
    drawnow;
end

if (isempty(param_file))
    error('ICA parameter file is not selected');
end

if (ischar(param_file))
    load(param_file);
    
    if ~exist('sesInfo', 'var')
        error('Not a valid parameter file');
    end
    outputDir = fileparts(param_file);
else
    sesInfo = param_file;
    outputDir = sesInfo.outputDir;
end

if (isempty(outputDir))
    outputDir = pwd;
end

try
    modalityType = sesInfo.modality;
catch
    modalityType = icatb_get_modality;
end

if (~strcmpi(modalityType, 'fmri'))
    warning('!!!Timecourse post-processing will be done only for fMRI modality');
    return;
end

writeInfo = 0;
try
    writeInfo = TIMECOURSE_POSTPROCESS.write;
catch
end

if (~writeInfo)
    return;
end

if (isfield(sesInfo, 'TR'))
    TR = sesInfo.TR;
else
    if (isempty(EXPERIMENTAL_TR))
        warning('!!!Experimental TR variable (EXPERIMENTAL_TR) in seconds is missing.');
        return;
    end
    TR = EXPERIMENTAL_TR;
end

if (length(TR) == 1)
    TR = repmat(TR, 1, sesInfo.numOfSub);
else
    if (length(TR) ~= sesInfo.numOfSub)
        error('Length of TR must match the number of subjects');
    end
end

% Spectra info
tapers = [3, 5];
sampling_frequency = 1/min(TR);
frequency_band = [0, 1/(2*min(TR))];

% FNC (despike timecourses and High freq cutoff in Hz)
despike_tc = 1;
cutoff_frequency = 0.15;

%% Write results
if (writeInfo)
    
    % SEPCTRA PARAMS (tapers, sampling_freq, frequency_band)
    try
        tapers = TIMECOURSE_POSTPROCESS.spectra.tapers;
    catch
    end
    
    %     try
    %         sampling_frequency = TIMECOURSE_POSTPROCESS.spectra.sampling_frequency;
    %     catch
    %     end
    %
    %     try
    %         frequency_band = TIMECOURSE_POSTPROCESS.spectra.frequency_band;
    %     catch
    %     end
    
    % FNC PARAMs (Despike timecourses and high frequency cutoff in Hz)
    try
        despike_tc = TIMECOURSE_POSTPROCESS.fnc.despike_tc;
    catch
    end
    
    try
        cutoff_frequency = TIMECOURSE_POSTPROCESS.fnc.cutoff_frequency;
    catch
    end
    
    outputFile = fullfile(outputDir, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
    
    %% Uncompress files
    subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', sesInfo.numOfSub, 'numOfSess', sesInfo.numOfSess, 'flagTimePoints', ...
        sesInfo.flagTimePoints);
    sesInfo.outputDir = outputDir;
    fileIn = dir(fullfile(outputDir, [sesInfo.calibrate_components_mat_file, '*.mat']));
    filesToDelete = {};
    if (length(fileIn) ~= sesInfo.numOfSub*sesInfo.numOfSess)
        disp('Uncompressing subject component files ...');
        filesToDelete = icatb_unZipSubjectMaps(sesInfo, subjectICAFiles);
    end
    
    %% Spectra
    spectra_params = struct('tapers', tapers, 'Fs', sampling_frequency, 'fpass', frequency_band);
    countS = 0;
    fprintf('\nComputing spectra and FNC correlations of all subjects and sessions components ...\n');
    if (despike_tc)
        disp('Timecourses will be despiked when computing FNC correlations...');
    end
    
    if (min(cutoff_frequency) > 0)
        disp(['Timecourses will be filtered when computing FNC correlations using HF cutoff of ', num2str(cutoff_frequency), ' Hz ...']);
    end
    
    minTpLength = min(sesInfo.diffTimePoints);
    if (~all(TR == min(TR)))
        tmpTR = TR;
        tmpTR = repmat(tmpTR(:)', sesInfo.numOfSess, 1);
        tmpTR = tmpTR(:)';
        ratiosTR = (tmpTR(:)')./min(tmpTR);
        [numN, denN] = rat(ratiosTR);
        chkTp = ceil((sesInfo.diffTimePoints(:)'.*numN)./denN);
        minTpLength = min(chkTp);
    end
    
    kurt_tc = zeros(sesInfo.numOfSub, sesInfo.numOfSess, sesInfo.numComp);
    kurt_ic = zeros(sesInfo.numOfSub, sesInfo.numOfSess, sesInfo.numComp);
    
    for nSub = 1:sesInfo.numOfSub
        for nSess = 1:sesInfo.numOfSess
            countS = countS + 1;
            timecourses = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'detrend_no', ...
                DETRENDNUMBER, 'subject_ica_files', subjectICAFiles);
            
            kvals = kurt(timecourses);
            
            kurt_tc(nSub, nSess, :) = kvals;
            
            
            % Interpolate timecourses if needed for variable TRs across
            % subjects
            if (~all(TR == min(TR)))
                interpFactor = TR(nSub)/min(TR);
                [num, denom] = rat(interpFactor);
                timecourses = resample(timecourses, num, denom);
            end
            
            %timecourses = timecourses(1:min(sesInfo.diffTimePoints), :);
            [temp_spectra, freq] = icatb_get_spectra(timecourses(1:minTpLength, :)', min(TR), spectra_params);
            temp_spectra = temp_spectra./repmat(sum(temp_spectra, 2), [1, size(temp_spectra, 2)]);
            temp_spectra = temp_spectra';
            if (countS == 1)
                spectra_tc_all = zeros(sesInfo.numOfSub, sesInfo.numOfSess, size(temp_spectra, 1), sesInfo.numComp);
                fnc_corrs_all = zeros(sesInfo.numOfSub, sesInfo.numOfSess, sesInfo.numComp, sesInfo.numComp);
            end
            spectra_tc_all(nSub, nSess, :, :) = temp_spectra;
            
            % despike
            if (despike_tc)
                timecourses = icatb_despike_tc(timecourses, min(TR));
            end
            
            % Filter
            if (min(cutoff_frequency) > 0)
                timecourses = icatb_filt_data(timecourses, min(TR), cutoff_frequency);
            end
            
            c = icatb_corr(timecourses);
            c(1:size(c, 1) + 1:end) = 0;
            c = icatb_r_to_z(c);
            fnc_corrs_all(nSub, nSess, :, :) = c;
        end
    end
    
    if (sesInfo.numOfSub*sesInfo.numOfSess == 1)
        spectra_tc_all = squeeze(spectra_tc_all);
        fnc_corrs_all = squeeze(fnc_corrs_all);
    end
    
    save(outputFile, 'spectra_tc_all', 'freq', 'fnc_corrs_all');
    
    %     %% FNC
    %     fprintf('\nComputing FNC correlations ...\n');
    %
    %     if (despike_tc)
    %         disp('Timecourses will be despiked ...');
    %     end
    %
    %     if (cutoff_frequency > 0)
    %         disp(['Timecourses will be filtered using HF cutoff of ', num2str(cutoff_frequency), ' Hz ...']);
    %     end
    
    %     countS = 0;
    %     fprintf('\n');
    %     for nSub = 1:sesInfo.numOfSub
    %         for nSess = 1:sesInfo.numOfSess
    %             countS = countS + 1;
    %             timecourses = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'tc', 'truncate_tp', 1, 'detrend_no', ...
    %                 DETRENDNUMBER, 'subject_ica_files', subjectICAFiles);
    %             for nC = 1:sesInfo.numComp
    %                 tmp = timecourses(:, nC);
    %                 if (despike_tc)
    %                     tmp = icatb_despike_tc(tmp, TR);
    %                 end
    %                 if (cutoff_frequency > 0)
    %                     tmp = icatb_filt_data(tmp, TR, cutoff_frequency);
    %                 end
    %                 timecourses(:, nC) = tmp;
    %             end
    %             if (countS == 1)
    %                 fnc_corrs_all = zeros(sesInfo.numOfSub, sesInfo.numOfSess, sesInfo.numComp, sesInfo.numComp);
    %             end
    %             c = icatb_corr(timecourses);
    %             c(1:size(c, 1) + 1:end) = 0;
    %             c = icatb_r_to_z(c);
    %             fnc_corrs_all(nSub, nSess, :, :) = c;
    %         end
    %     end
    
    %     if (sesInfo.numOfSub*sesInfo.numOfSess == 1)
    %         fnc_corrs_all = squeeze(fnc_corrs_all);
    %     end
    %
    %     save(outputFile, 'fnc_corrs_all', '-append');
    
    
    % mutual information between components
    countS = 0;
    fprintf('\n');
    for nSub = 1:sesInfo.numOfSub
        for nSess = 1:sesInfo.numOfSess
            countS = countS + 1;
            ic = icatb_loadComp(sesInfo, (1:sesInfo.numComp), 'subjects', nSub, 'sessions', nSess, 'vars_to_load', 'ic', 'subject_ica_files', subjectICAFiles);
            kvals = kurt(ic);
            kurt_ic(nSub, nSess, :) = kvals;
            tmp = icatb_compute_mi(ic');
            if (countS == 1)
                spatial_maps_MI = zeros(sesInfo.numOfSub, sesInfo.numOfSess, sesInfo.numComp, sesInfo.numComp);
            end
            spatial_maps_MI(nSub, nSess, :, :) = tmp;
            clear tmp;
        end
    end
    
    kurt_comp.ic = kurt_ic;
    kurt_comp.tc = kurt_tc;
    
    save(outputFile, 'spatial_maps_MI', 'kurt_comp', '-append');
    
    fprintf('Done\n\n');
    disp(['File ', outputFile, ' contains spectra and FNC correlations']);
    disp('1. spectra_tc_all - Timecourses spectra. spectra_tc_all variable is of dimensions subjects x sessions x spectral length x components');
    disp('2. fnc_corrs_all - FNC correlations transformed to fisher z-scores. fnc_corrs_all variable is of dimensions subjects x sessions x components x components');
    disp('3. spatial_maps_MI - Mutual information is computed between components spatially. spatial_maps_MI variable is of dimensions subjects x sessions x components x components');
    disp('4. kurt_comp - Kurtosis is computed on the spatial maps and timecourses. kurt_comp variable is of dimensions subjects x sessions x components');
    
    
    if (exist('filesToDelete', 'var') && ~isempty(filesToDelete))
        icatb_cleanupFiles(filesToDelete, outputDir);
    end
    
end



function k = kurt(x)
% Kurtosis

x = icatb_remove_mean(x);
s2 = mean(x.^2);
m4 = mean(x.^4);
k = (m4 ./ (s2.^2));