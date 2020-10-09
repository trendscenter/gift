function [V, Lambda, V_orig, Lambda_orig] = icatb_v_pca(data, firstEig, lastEig, saveEig, doTranspose, removeMean)
%% compute covariance matrix by removing the mean and then calculate the
% eigen values and eigen vectors

% Input: 1. data - 2D array, 1 dimension - voxels, 2 dimension - time points
% 2., 3. and 4.,  firstEig, lastEig and saveEig - scalars
% 5. doTranspose - 'transpose' - computes transpose

firstEig_t2 = firstEig; lastEig_t2 = lastEig; saveEig_t2=saveEig;

if (~exist('saveEig','var'))
    saveEig = 0;
end;

if ~exist('removeMean', 'var')
    removeMean = 'yes';
end

if (saveEig == 3),
    disp('Loading Previous EigenDecomposition 1');edit
    load eig_stuff1;
    firstEig = firstEig_t2;lastEig=lastEig_t2;saveEig=saveEig_t2;
    originalDimension = size(data, 1);
else
    originalDimension = size(data, 1);

    if exist('doTranspose', 'var')
        if strcmp(doTranspose, 'transpose')
            data = data';
        end
        originalDimension = size(data, 2);
    end

    disp(['Calculating Covariance: ' num2str(originalDimension) '^2']);

    %% Skip detrending if the removing the mean is already done
    rmMean = 0;
    if strcmpi(removeMean, 'yes')
        disp('Removing mean from the data ...');
        rmMean = 1;
    end

    %% Calculate covariance matrix
    covarianceMatrix = icatb_cov(data, rmMean);

    clear data;

    %% Calculate the eigenvalues and eigenvectors of covariance matrix.
    disp('Calculating eigendecomposition');

    condCovMat = cond(covarianceMatrix);

    disp(['Condition number of covariance matrix is ', num2str(condCovMat)]);

    [V, Lambda] = eig(covarianceMatrix, 'nobalance');

    %% Sort the eigenvalues - decending.
    disp('Sorting eigenvalues');
    [eigenvalues ind] = sort(diag(Lambda));
    eigenvalues = flipud(eigenvalues);
    ind = flipud(ind);

    if (saveEig == 4),
        disp('Saving Eigen Decomposition 1');
        save eig_stuff1 Lambda V covarianceMatrix eigenvalues;
        clear covarianceMatrix;
    end

end

%estimate the rank (Too long for tica!)...
%maxLastEig = rank(covarianceMatrix);
maxLastEig = lastEig;

disp('Selecting Desired Eigenvalues');
% Remove the smaller eigenvalues
%disp('Removing smaller eigenvalues');
if lastEig < originalDimension
    lowerLimitValue = eigenvalues(lastEig);
else
    lowerLimitValue = eigenvalues(originalDimension);
end
lowerColumns = (diag(Lambda) >= lowerLimitValue);

%disp('Removing larger eigenvalues');
% Remove the larger eigenvalues
higherLimitValue = eigenvalues(firstEig);
higherColumns = diag(Lambda) <= higherLimitValue;

% Combine the results from above
selectedColumns = lowerColumns & higherColumns;

% Select the colums which correspond to the desired range
% of eigenvalues.

%% Fixed this statement in the updates of GroupICATv2.0b (Mentioned by Jean
% liu)
%V_orig = V(ind,ind);

V_orig = V(:, ind);

V = icatb_v_selcol(V,selectedColumns);
Lambda_orig = eigenvalues;
Lambda = icatb_v_selcol(icatb_v_selcol(Lambda, selectedColumns)',selectedColumns);

sumAll=sum(eigenvalues);
sumUsed=sum(diag(Lambda));
retained = (sumUsed/sumAll)*100;
fprintf('%g%% of (non-zero) eigenvalues retained.\n', retained);
