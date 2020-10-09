function icatb_spm_avg_runs(parameter_file)
% Average runs for subjects for all components using spm_imcalc

% Load defaults
icatb_defaults;
global PARAMETER_INFO_MAT_FILE;
global MEAN_AN3_FILE;
global COMPONENT_NAMING;
global SPM_STATS_AVG_RUNS;

if (SPM_STATS_AVG_RUNS)
    
    filterP = ['*', PARAMETER_INFO_MAT_FILE, '*.mat'];
    if ~exist('parameter_file', 'var') || isempty(parameter_file)
        parameter_file = icatb_selectEntry('typeEntity', 'file', 'typeSelection', 'single', 'filter', ...
            filterP, 'title', 'Select a valid parameter file');
    end
    
    if isempty(parameter_file)
        error('Parameter file is not selected');
    end
    
    load(parameter_file);
    
    if ~exist('sesInfo', 'var')
        error(['Selected file: ', parameter_file, ' is not a valid parameter file']);
    end
    
    outputDir = fileparts(parameter_file);
    
    % Get necessary fields
    numOfSub = sesInfo.numOfSub;
    numOfSess = sesInfo.numOfSess;
    numComp = sesInfo.numComp;
    zipContents = sesInfo.zipContents;
    prefix = sesInfo.userInput.prefix;
    
    % Check if number of sessions is greater than 1
    if numOfSess > 1
        
        spmPath = which('spm.m');
        
        if isempty(spmPath);
            error('SPM does not exist on MATLAB path');
        end
        
        verNum = str2num(strrep(lower(spm('ver')), 'spm', ''));
        
        if (verNum < 5)
            error('This utility works with SPM5 and higher');
        end
        
        % Subject ICA files
        subjectICAFiles = icatb_parseOutputFiles('icaOutputFiles', sesInfo.icaOutputFiles, 'numOfSub', numOfSub, 'numOfSess', numOfSess, ...
            'flagTimePoints', sesInfo.flagTimePoints);
        
        disp('Checking subject component maps ...');
        filesToDelete = icatb_unZipSubjectMaps(sesInfo, subjectICAFiles);
        % End for checking if zip files are present or not
        disp('Done checking subject component maps');
        
        disp('Reading headers ...');
        % Loop over subjects
        for nSub = 1:numOfSub
            % Loop over sessions
            for nSess = 1:numOfSess
                compFiles = icatb_rename_4d_file(subjectICAFiles(nSub).ses(nSess).name);
                subjectICAFiles(nSub).ses(nSess).name = compFiles;
                clear compFiles;
            end
            % End loop over sessions
        end
        % End loop over subjects
        disp('Done reading headers');
        
        % Image calculation formula
        imCalcFormula = imcalc_formula(numOfSess);
        
        meanOverRunsDir = [prefix, '_avg_runs'];
        if ~exist(fullfile(outputDir, meanOverRunsDir), 'dir')
            mkdir(outputDir, meanOverRunsDir);
        end
        
        % Loop over component numbers
        outFiles = cell(numComp, numOfSub);
        % Loop over components
        for compNum = 1:numComp
            % Loop over subjects
            for nSub = 1:numOfSub
                files = cell(numOfSess, 1);
                % Loop over sessions
                for nSess = 1:numOfSess
                    files{nSess} = deblank(subjectICAFiles(nSub).ses(nSess).name(compNum, :));
                end
                % End loop over sessions
                files = strcat(outputDir, filesep, str2mat(files));
                out_file_name = [prefix, MEAN_AN3_FILE, '_sub_', icatb_returnFileIndex(nSub), COMPONENT_NAMING, icatb_returnFileIndex(compNum), '.img'];
                disp(['Calculating mean over runs for subject ', icatb_returnFileIndex(nSub), ' component ', icatb_returnFileIndex(compNum)]);
                fprintf('\n');
                % Compute mean image
                spm_imcalc_ui(files, fullfile(outputDir, meanOverRunsDir, out_file_name), imCalcFormula, {0, 0, 16, 0});
                % Store the file names
                outFiles{compNum, nSub} = out_file_name;
            end
            % End loop over subjects
        end
        % End loop over components
        clear subjectICAFiles;
        
        % Store files field in sesInfo
        sesInfo.spm_stats.avg_runs_files = outFiles;
        sesInfo.spm_stats.avg_runs_dir = meanOverRunsDir;
        clear outFiles;
        
        % Save the information in a parameter file
        icatb_save(parameter_file, 'sesInfo');
        
        icatb_cleanupFiles(filesToDelete, outputDir);
        
    end
    % End for checking if the number of sessions is greater than 1
    
end

function imCalcFormula = imcalc_formula(numOfSess)
% imcalc formula

imCalcFormula = '(';
% Loop over sessions
for nSess = 1:numOfSess
    if (nSess == 1)
        imCalcFormula = [imCalcFormula, 'i', num2str(nSess)];
    else
        imCalcFormula = [imCalcFormula, '+i', num2str(nSess)];
    end
end
% End loop over sessions

if (nSess == numOfSess)
    imCalcFormula = [imCalcFormula, ')/', num2str(numOfSess)];
end