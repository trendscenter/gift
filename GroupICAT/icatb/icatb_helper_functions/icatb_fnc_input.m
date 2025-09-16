function [fnc_matrix] = icatb_fnc_input(inputData)
%% Return FNC matrix
%

fnc_variable_mat_file = inputData.fnc_variable_mat_file;


%% Read data
input_data_file_patterns = cellstr(inputData.input_data_file_patterns);
if (numel(input_data_file_patterns) > 1)
    countS = 0;
    for nSub = 1:size(input_data_file_patterns, 1)
        for nSess = 1:size(input_data_file_patterns, 2)
            countS = countS + 1;
            current_file_name = input_data_file_patterns{nSub, nSess};
            [pathstr, fileN, extn] = fileparts(current_file_name);
            fileList = icatb_listFiles_inDir(pathstr, [fileN, extn]);
            fileList = icatb_fullFile('files', fileList, 'directory', pathstr);
            fileList = deblank(fileList(1, :));
            fnc_results = load(fileList);
            if (strcmpi(extn, '.mat'))
                fnc_mat_files = fnc_results.(fnc_variable_mat_file);
            else
                fnc_mat_files= fnc_results;
            end
            if (countS == 1)
                fnc_matrix = zeros(size(input_data_file_patterns, 1), size(input_data_file_patterns, 2), size(fnc_mat_files, 1), size(fnc_mat_files, 2));
            end
            fnc_matrix(nSub, nSess, :, :) = fnc_mat_files;
        end
    end
else
    % check if parameter file
    [pathstr, fileN, extn] = fileparts(input_data_file_patterns{1});
    fileList = icatb_listFiles_inDir(pathstr, [fileN, extn]);
    fileList = icatb_fullFile('files', fileList, 'directory', pathstr);
    if (size(fileList, 1) == 1)
        % Check if it is a parameter file
        load(fileList);
        if (exist('sesInfo', 'var'))
            postprocessFile = fullfile(pathstr, [sesInfo.userInput.prefix, '_postprocess_results.mat']);
            load(postprocessFile);
            if exist('fnc_corrs_all', 'var')
                fnc_matrix = fnc_corrs_all;
            else
                postProcessFiles = icatb_fullFile('directory', pathstr, 'files', char(postProcessFiles));
                fnc_mat_files = cell(1, size(postProcessFiles, 1));
                for nFile = 1:size(postProcessFiles, 1)
                    c_file_name = deblank(postProcessFiles(nFile, :));
                    load(c_file_name);
                    fnc_mat_files{nFile} = fnc_corrs;
                    clear fnc_corrs;
                end
                fnc_matrix = permute(cat(4, fnc_mat_files{:}), [4, 1, 2, 3]);
            end
        else
            % Load FNC mat file
            fnc_results = load(fileList, fnc_variable_mat_file);
            fnc_matrix = fnc_results.(fnc_variable_mat_file);
        end
        
    else
        % Load fnc files
        fnc_mat_files = cell(1, size(fileList, 1));
        for nFile = 1:size(fileList, 1)
            c_file_name = deblank(fileList(nFile, :));
            [~, fileN, extn] = fileparts(input_data_file_patterns{1});
            fnc_results = load(c_file_name);
            if (strcmpi(extn, '.mat'))
                fnc_mat_files{nFile} = fnc_results.(fnc_variable_mat_file);
            else
                fnc_mat_files{nFile} = fnc_results;
            end
        end
        fnc_matrix = permute(cat(3, fnc_mat_files{:}), [3, 1, 2]);
        fnc_matrix = reshape(fnc_matrix, [size(fnc_matrix, 1), 1, size(fnc_matrix, 2), size(fnc_matrix, 3)]);
        
    end
    
end

contrast_vector = [];
try
    contrast_vector = inputData.contrast_vector;
catch
end

if (~isnumeric(contrast_vector))
    contrast_vector = str2num(contrast_vector);
end

if ~isempty(contrast_vector)
    if length(contrast_vector) ~= size(fnc_matrix, 2)
        error('contrast vector doesn''t match the sessions');
    end
end


% contrast vector
if isempty(contrast_vector)
    fnc_matrix = squeeze(mean(fnc_matrix, 2));
else
    
    fnc_agg = zeros(size(fnc_matrix, 1), size(fnc_matrix, 3), size(fnc_matrix, 4));
    for nSess = 1:size(fnc_matrix, 2)
        fnc_agg = fnc_agg + (contrast_vector(nSess)*squeeze(fnc_matrix(:, nSess, :, :)));
    end
    
    fnc_matrix = fnc_agg;
    
end

clear sesInfo;