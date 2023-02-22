function [icasig, A, betas] = icatb_percentSignalChangeN(origData, icasig, A, mask_ind)
% scales to percent signal change

A = A';

if ~exist('mask_ind', 'var')
    mask_ind = find(icasig(1, :) ~= 0);
end

% Loop over components
[icasig] = icatb_convertImageToZScores(icasig);

icasigN = icasig';

if size(icasig, 2) ~= length( mask_ind)
    icasigN = icasig(:, mask_ind)';
end

if size(origData, 1) ~= length( mask_ind)
    origData = origData(mask_ind, :);
end


disp('Using a new method to calculate percent signal change');

% Initialise betas
betas = zeros(size(A, 1), [size(icasigN, 2) + 1]);

numComp = size(A, 2);

% Loop over time points
for ii = 1:size(A, 1)
   % fit component data using raw data as observation
    
    obs = origData(:, ii); % observation of dimensions voxels by 1  
    X = [icasigN, ones(size(icasigN, 1), 1)]; % form model matrix
    % 
    [betaWeights, R2] = icatb_regress(obs, X);    
    
    betas(ii, :) = betaWeights;
    
end

A = betas(:, 1:numComp);

% % Loop over components
% for nn = 1:size(A, 2)
%     A(:, nn) = A(:, nn).*abs(betas(:, nn));
% end


%save('percent_signal_change.mat', 'betas', 'A');

%disp('Saved results of percent signal change');
