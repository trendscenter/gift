% noisecloud_classify - for use with extracted features to build a
% classifier with logistic regression, using the elastic net penalty.
% Input should be normalized features produced with noisecloud.m.

% REQUIREMENTS: glmnet package for matlab (included)

function RESULT = noisecloud_classify(features_norm,feature_labels, labels, varargin)


iterations = 1;
K = 10;

for i = 1:2:length(varargin)
    if (strcmpi(varargin{i}, 'iterations'))
        iterations = varargin{i + 1};
    elseif (strcmpi(varargin{i}, 'cross_validation'))
        K = varargin{i + 1};
    end
end

% Add glmnet package
addpath(genpath('3rdparty/glmnet_matlab'));

%     % Step 1: Ask user for network labels
%     [label_file,~] = spm_select(1,'mat','Select component labels (noise/real) .mat file ','','','.mat');
%     labels = load(label_file);  labels = labels.output;
%
%     % The labels are 1==noise, 0==not.  Change 0s to 2s
labels.decision(labels.decision == 0) = 2;

% Show imagesc of features and components (subnetworks)
figure;
imagesc(features_norm); colorbar
title([ num2str(size(features_norm,2)) ' Normalized Spatial and Temporal Features for ' num2str(size(features_norm,1)) ' Components ' ]);
ylabel('Components');
xlabel('Features');

% X: data; y: class labels; K = number of cross validation folds
% (usuallly K = 10); Lambda = [0.0625 0.125 0.25 0.5 1 2 4 8 16];
% type = 'binomial' or 'multinomial'

% Ask the user to input running parameteres
%iterations = input('Enter iterations (number of permutations) (default==1):');
%K = input('Enter K value for cross validation (default==10):');
type = 'binomial';
Lambda = [0.0625 0.125 0.25 0.5 1 2 4 8 16];
alphav = 1:-0.01:0.1;

fprintf('\n%s\n','Running GLMNET (logistic regression with elastic net');
fprintf('%s\n',[ 'K: ' num2str(K) ]);
fprintf('%s\n',[ 'Grid search for alpha between 0,1 and 1'   ]);
fprintf('%s\n\n',[ 'Type: ' type   ]);

RESULT = noisecloud_runGlmnetCV(features_norm,labels.decision,K,iterations,alphav,Lambda,type,feature_labels);

end