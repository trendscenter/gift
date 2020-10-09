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