function [class_labels, fit_mdl, result_nc_classifier] = noisecloud_run(training_opts, testing_opts, varargin)
%% Run noise cloud (requires SPM)
%
% Inputs:
% 1. training_opts: data structure
% a. TR - TR in seconds
% b. class_labels - Class labels of length equal to the number of
% components in training data. It has two levels where 1 corresponds to
% noise and 0 corresponds to network.
% c. sm - Spatial map file names
% d. tc - Timecourse file names (Enter one file name per subject)
% 2. testing_opts: data structure with the following fields:
%   a. TR - TR in seconds
%   b. sm - Spatial map file names
%   c. tc - Timecourse file names (Enter one file name per subject)
%
% varargin - Variable no of arguments
%
% Outputs:
%
% class_labels - Class labels of the testing data-set. 1 for noise and 0
% for network detected.
% fit_mdl - Model fit using the best lambda from Noise cloud classifier
% result_nc_classifier - Results from noise cloud classifier
%
% nc_class_labels.txt file containing labels (Noise or network) is written in the output directory specified.


outDir = pwd;
coregister_im = 1;
num_iterations = 1;
num_cross_validation = 10;
convert_to_z = 'yes';
for nV = 1:length(varargin)
    if (strcmpi(varargin{nV}, 'convert_to_z'))
        convert_to_z = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'outDir'))
        outDir = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'coregister'))
        coregister_im = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'iterations'))
        num_iterations = varargin{nV + 1};
    elseif (strcmpi(varargin{nV}, 'cross_validation'))
        num_cross_validation = varargin{nV + 1};
    end
end


file_names_testing = [];
file_names_training = [];

if (isfield(training_opts, 'sm') && ~isempty(training_opts.sm))
    
    % Form full file names of training data
    [file_names_training, dims1] = getFilenames(training_opts.sm);
    
    % Form full file names of testing data
    [file_names_testing, dims2] = getFilenames(testing_opts.sm);
    
    % Coregister testing data to training data if the dimensions are different
    if ~all(dims1 == dims2)
        otherImages = '';
        if (length(file_names_testing) > 1)
            otherImages = char(file_names_testing(2:end));
        end
        disp(['Coregistering ', file_names_testing{1}, ' to functional image ', file_names_training{1}]);
        noisecloud_spm_coregister(file_names_training{1}, file_names_testing{1}, otherImages, outDir);
        newFileNames = cell(1, length(file_names_testing));
        for nF = 1:length(file_names_testing)
            cF = file_names_testing{nF};
            [pp, fN, extn] = fileparts(cF);
            newFileNames{nF} = fullfile(outDir, ['r', fN, extn]);
        end
        file_names_testing = newFileNames;
    end
    
end

training_tc = [];
testing_tc = [];
regress_training_cov = [];
regress_testing_cov = [];
if (isfield(training_opts, 'tc') && ~isempty(training_opts.tc))
    training_tc = training_opts.tc;
    testing_tc = testing_opts.tc;
    try
        regress_training_cov = training_opts.regress_cov;
    catch
    end
    try
        regress_testing_cov = testing_opts.regress_cov;
    catch
    end
end

orig_data = '';
try
    orig_data = testing_opts.orig_data;
    orig_data(1);
    orig_data = cellstr(orig_data);
catch
end

[training_tc, testing_tc] = doRegressTimeCourses(training_tc, regress_training_cov, testing_tc, regress_testing_cov, outDir);

%% Get the features of training data
disp('Computing spatial and temporal features of training data ...');
[features_norm, feature_labels] = noisecloud(training_opts.TR, file_names_training, training_tc, 'convert_to_z', convert_to_z, 'outDir', outDir, 'coregister', coregister_im);

if (isnumeric(training_opts.class_labels))
    tmpLa = training_opts.class_labels;
else
    tmpLa = load(training_opts.class_labels);
end

labels.decision = tmpLa(:);
disp('Building classifier using logistic regression ...');
result_nc_classifier = noisecloud_classify(features_norm, feature_labels, labels, 'iterations', num_iterations, 'cross_validation', num_cross_validation);

%% Use the best lambda parameter to train the data again
options = glmnetSet;
options.alpha = result_nc_classifier.opt_alpha;
options.lambda = result_nc_classifier.opt_lambda;
labels.decision(labels.decision == 0) = 2;
fit_mdl = glmnet(features_norm, labels.decision, 'binomial', options);

%% Use testing data-set to predict the labels
disp('Computing spatial and temporal features of testing data ...');
[features_norm_test, feature_labels_test] = noisecloud(testing_opts.TR, file_names_testing, testing_tc, 'convert_to_z', convert_to_z, 'outDir', outDir, ...
    'coregister', 0);
disp('Predicting ...');
class_labels = glmnetPredict(fit_mdl, 'class', features_norm_test);
class_labels(class_labels == 2) = 0;
fprintf('\nDone\n');

%% Write info in text file
flags_testing = repmat({'Noise'}, length(file_names_testing), 1);
flags_testing(class_labels == 0) = {'Network'};
output_file_name = fullfile(outDir, 'nc_class_labels.txt');
cl_to_write = [char(file_names_testing), repmat(' ', length(file_names_testing), 1), char(flags_testing)];
dlmwrite(output_file_name, cl_to_write, '');

%% Optional remove noise from the testing data-sets
if (~isempty(orig_data))
    if (~isempty(file_names_testing))
        removeNoiseData(orig_data, testing_opts.sm, class_labels, outDir);
    end
end



function [file_names, dims] = getFilenames(network_paths)
% Form full filenames from 4D nifti data
%
VV = spm_vol(network_paths);
dims = VV(1).dim(1:3);
file_names = cell(length(VV), 1);
for nV = 1:length(VV)
    file_names{nV} = [VV(nV).fname, ',', num2str(VV(nV).n(1))];
end



function [training_tc, testing_tc] = doRegressTimeCourses(training_tc, regress_training_cov, testing_tc, regress_testing_cov, outDir)
% Regress timecourses

if (~isempty(training_tc))
    training_tc = cellstr(training_tc);
    if (~isempty(regress_training_cov))
        regress_training_cov = cellstr(regress_training_cov);
        training_tc = regressTC(training_tc, regress_training_cov, outDir);
    end
    training_tc = char(training_tc);
end


if (~isempty(testing_tc))
    testing_tc = cellstr(testing_tc);
    if (~isempty(regress_testing_cov))
        regress_testing_cov = cellstr(regress_testing_cov);
        testing_tc = regressTC(testing_tc, regress_testing_cov, outDir);
    end
    testing_tc = char(testing_tc);
end



function training_tc = regressTC(training_tc, regress_training_cov, outDir)

for nS = 1:length(training_tc)
    cF = training_tc{nS};
    [p, fn, extn] = fileparts(cF);
    commaPos = find(extn == ',');
    if (~isempty(commaPos))
        extn = extn(1:commaPos(1)-1);
    end
    outName = fullfile(outDir, ['r_', fn, extn]);
    
    if (~strcmpi(extn, '.nii') && ~strcmpi(extn, '.img'))
        % load as ascii
        tc = icatb_load_ascii_or_mat(cF);
    else
        % load as nifti
        V = spm_vol(cF);
        tc = spm_read_vols(V);
    end
    covTc = icatb_load_ascii_or_mat(regress_training_cov{nS});
    if (numel(covTc) == length(covTc))
        covTc = covTc(:);
    end
    tc = detrend(tc, 0);
    betas = pinv(covTc)*tc;
    tc = tc - covTc*betas;
    if (~strcmpi(extn, '.nii') && ~strcmpi(extn, '.img'))
        % save as ascii file
        save(outName, 'tc', '-ascii');
    else
        % save as nifti file
        V.fname = outName;
        spm_write_vol(V, tc);
    end
    training_tc{nS} = outName;
end


function removeNoiseData(orig_data, file_names_testing, class_labels, outDir)
% Remove noise components from the data

file_names_testing = cellstr(file_names_testing);

endT = 0;
for nO = 1:length(orig_data)
    currentFile = orig_data{nO};
    [p, fn, extn] = fileparts(currentFile);
    tmp_files = dir(currentFile);
    if (isempty(tmp_files))
        error('File %s names doesn''t exist with the pattern or file name', currentFile);
    end
    p = icatb_fullFile('directory', p, 'files', char(tmp_files.name));
    files = icatb_rename_4d_file(p);
    V = icatb_spm_vol(deblank(files));
    Vcomp = icatb_spm_vol(file_names_testing{nO});
    mask = find(abs(icatb_spm_read_vols(Vcomp(1))) > eps);
    data = icatb_read_data(files, [], (1:prod(V(1).dim(1:3))));
    startT = endT + 1;
    endT = endT + length(Vcomp);
    comp_flags = class_labels(startT:endT);
    noise_comps = find(comp_flags == 1);
    if (isempty(noise_comps))
        disp(['No noise components in ', file_names_testing{nO}, '. Not removing components']);
        continue;
    end
    sm = icatb_remove_mean(icatb_read_data(file_names_testing{nO}, [], mask));
    dataN = data(mask, :);
    meanData = mean(dataN);
    % remove the mean
    dataN = icatb_remove_mean(dataN, 0);
    A = (pinv(sm)*dataN)';
    S = sm';
    A_R = A; S_R = S;
    A_R(:, noise_comps) = []; S_R(noise_comps, :) = [];
    % compute the new data (contains the removed articfactual)
    disp('Removing the noise components from the data ...');
    dataN = A_R*S_R;
    clear  A S A_R S_R;
    disp('Adding mean back to the data ...');
    dataN = dataN';
    dataN = bsxfun(@plus, dataN, meanData);
    %dataN = zeros(prod(V(1).dim(1:3)), size(data, 2));
    data(mask, :) = dataN;
    clear dataN;
    data = reshape(data, [V(1).dim(1:3), size(data, 2)]);
    newDir = fullfile(outDir, 'noise_cloud_cleaned_data');
    if (exist(newDir, 'dir') ~= 7)
        mkdir(newDir);
    end
    [~, tmpP, extn] = fileparts(V(1).fname);
    outFile = fullfile(newDir, ['R_', tmpP, '.nii']);
    icatb_write_nifti_data(outFile, V, data, 'Cleaned data');
    clear data;
end

fprintf('Done\n');
