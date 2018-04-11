% Output directory to place the coregistered images 
outputDir = 'C:\Users\srachakonda\Desktop\test_nc';

%% labels matrix - 1 for noise and 0 for networks
comp_network_names = {'BG', 21;% Basal ganglia 21st component
'AUD', 17;                   % Auditory 17th component
'SM', [7 23 24 29 38 56];    % Sensorimotor comps
'VIS', [39 46 48 59 64 67];  % Visual comps
'DMN', [25 50 53 68];        % DMN comps
'ATTN', [34 52 55 60 71 72]; % ATTN Comps
'FRONT', [20 42 47 49]};     % Frontal comps
vals = comp_network_names(:,2);
vals = [vals{:}];
labels = ones(75, 1);
labels(vals) = 0;


%% Training data 
training_opts.sm = 'C:\Users\srachakonda\Desktop\test_nc\dat\rest_hcp_mean_component_ica_s_all_.nii'; % Spatial maps
training_opts.tc = 'C:\Users\srachakonda\Desktop\test_nc\dat\rest_hcp_mean_timecourses_ica_s_all_.nii'; % Timecourses

% Specify any covariates that need to be regressed out from the data in a
% cell array. 
training_opts.regress_cov = [];

% TR of training data
training_opts.TR = 2;
training_opts.class_labels = labels; % Labels 

%% Testing data
testing_opts.sm = 'E:\test_GIFT\new_version\Multiple_sub_Multiple_sess\Visuo_sub001_component_ica_s1_.nii';
testing_opts.tc = 'E:\test_GIFT\new_version\Multiple_sub_Multiple_sess\Visuo_sub001_timecourses_ica_s1_.nii';

% Specify any covariates that need to be regressed out from the data in a
% cell array. 
testing_opts.regress_cov = [];

% TR of testing data 
testing_opts.TR = 1;

[class_labels, fit_mdl, result_nc_classifier] = noisecloud_run(training_opts, testing_opts, 'convert_to_z', 'yes', 'outDir', outputDir, 'coregister', 0, ...
    'iterations', 1, 'cross_validation', 10);