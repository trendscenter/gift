function [L,randmatInfo] = icatb_dist2RandMat(inpMat,matType,nIter,matParams)
%
% Created by Victor M. Vergara
%
% First Version: 2016-11-08
% 2018-01-07: In this January version we updated the normalization of the
% metric L. Originally it was L = ((Z/svcov)*Z')./rank_of_M;, but this was
% not working for rank deficient matrices. In this version we changed to
% the smallest dimmension of the matrix L = ((Z/svcov)*Z')./small_dim_of_M;
%
% Distance to Random Matrices
%
% [L,randmatInfo] = icatb_dist2RandMat(inpMat,matType,nIter,matParams)
%
% Finds the distance of matrix inpMat to a random matrix with the same
% properties as inpMat including size, matrix element mean and matrix
% element variance.
%
% Input:
%
%	inpMat : the matrix to be tested
%
%	matType: can assume the values 0, 1 or 2
%		The code doesn't detect symmetry, so it needs to be indicated by
%		the input variable matType.
% matType=0 indicate M is not symmetric and can be square (Default)
% matType=1 indicate M is symmetric (all main diagonal elements are random)
% matType=2 indicate M is a correlation matrix(main diagonal elements = 1)
%
%	  nIter: number of iterations to estimate means and covariance of SVs.
%		     if matParams is defined, then nIter is not used
%			(Default:	nIter = 1000)
%
% matParams: cell must contain mean SVs and covariance SVs
%		matParams{1} = mean singular values previously calculated
%		matParams{2} = covariance matrix of the singular values
%
%           if matParams is not defined then the code will estimate the
%           parameters, but will take longer depending on nIter.
%
%% Since 2018-01-07 we always normalize by the smallest dimension
[n,m] = size(inpMat);
small_dim_of_inpMat = min([n,m]);
%% Make sure matType is defined
if ~exist('matType','var')
    matType = 0;
end
if (matType<0) || (matType>2)
    warning(['Invalid value for matType=',num2str(matType),	...
        ', but setting matType=0 [general matrix case]']);
    matType = 0;
end
%% Get the covariance matrix and the mean of all SVs
if exist('matParams','var')
    svmean = matParams{1};
    svcov = matParams{2};
    elemean= mean(inpMat(:));
    elemstd= std(inpMat(:));
else
    %% Make sure the number of iterations is defined
    if ~exist('nIter','var')
        nIter = 100000;
        warning(['Number of Iterations Set to ',num2str(nIter)]);
    end
    [svmean,svcov,elemean,elemstd] = ...
        calculate_rand_mat_stats([n,m],inpMat,matType,nIter);
end
%% Calculate Randomness Metric L
sv = svd(inpMat)';
Z = (sv - svmean);
L = ((Z/svcov)*Z')./small_dim_of_inpMat;
%% Calculate the p-value
pval = chi2cdf( L*small_dim_of_inpMat,	small_dim_of_inpMat, 'upper');
%% Fill out output information structure
randmatInfo.L = L;
randmatInfo.pval = pval;
randmatInfo.singularValuesMean = svmean;
randmatInfo.singularValuesCovariance = svcov;
randmatInfo.matrixElementsMean = elemean;
randmatInfo.matrixElementsVariance = elemstd;

%% This one is needed to calculate the statistics of random matrices
function [svmean,svcov,elemean,elemstd] = calculate_rand_mat_stats(sZ,M,matType,Niter)
%% Attend to each case according to variable matType
% Elements mean and std must be calculated differently for each type
switch matType
    case 0
        V = M(:);
    case 1
        V = M(:);
    case 2
        V = M;for kk=1:size(M,1);V(kk,kk)=1;end;
        V = V(:);
    otherwise
        error(['Invalid value of s=',num2str(s),', but setting s=0']);
end
% Estimate the random matrices
elemean = mean(V);
elemstd =  std(V);
[svmean,svcov] = vvrandsvstat(sZ,elemean,elemstd, Niter, matType );



function [svmean,svcov,svs] = vvrandsvstat(sZ, xmean, xstd,  Niter, matype)
%
% Created by Victor M. Vergara
%
% First Version: 2016-11-08
%
% sv = randsvspec( N, M, xmean, xstd, Niter )
%
% Singular Values Statistics for random matrices
%
% Random Matrix of size N x M
% The elements of the matrix are Gaussian with (xmean, xstd)
%
% Niter: number of iterations to estimate
%
% matype: type of matrix
% The code doesn't detect symmetry, so it needs to be indicated by the
% input variable s.
% s=0 indicate M is not symmetry (Default)
% s=1 indicate M is symmetric(all main diagonal elements are random)
% s=2 indicate M is a correlation matrix(all main diagonal elements are 1)
%
%
%% Complete missing values
if ~exist('matype','var')
    matype = 0;
end
if ~exist('Niter','var')
    Niter = 100000;
end
if ~exist('xmean','var')
    xmean = 0;
end
if ~exist('xstd','var')
    xstd = 1.0;
end
%% process each matrix type acording to the case
switch matype
    case 0
        svs = nonsymrandsv(sZ(1),sZ(2), xmean, xstd, Niter );
    case 1
        svs = symrandsv(sZ(1), xmean, xstd, Niter);
    case 2
        svs = corsymrandsv(sZ(1), xmean, xstd,Niter);
    otherwise
        error(['vvrandsvstat:Invalid value of s=',num2str(s), ...
            ', but setting s=0']);
end
svmean = mean(svs);
svcov  = cov(svs);

%% Create a symmetric matrix by copying the upper diagnoal into the
% lower diagonal and ignores the lower diagonal elements
function H = dothesym(A)
H = triu(A); % Set lower triangular elements to 0
T = triu(A,1)'; % Transpose T
for kk=1:size(T,1);T(kk,kk)=0;end;% set main diag to 0
H = H + T; % copy elements above he diagonal to the lower

%% Set the main diagnoal to 1s
function H = setdiag2one(A)
H = A;
for kk=1:size(A,1)
    H(kk,kk) = 1;
end
%% Estimate sv_mean and sv_std for the symmetric correlation case
% the main diagonal is all ones
function svs = corsymrandsv(Nsv, xmean, xstd, Niter)
svs = zeros(Niter,Nsv);
for iter=1:Niter
    H = xstd.*randn(Nsv) + xmean;
    H = dothesym(H);
    H = setdiag2one(H);
    svs(iter,:) = svd(H);
end

%% Estimate sv_mean and sv_std for the symmetric case
function svs = symrandsv(Nsv, xmean, xstd, Niter)
svs = zeros(Niter,Nsv);
for iter=1:Niter
    H = xstd.*randn(Nsv) + xmean;
    H = dothesym(H);
    svs(iter,:) = svd(H);
end

%% Estimate sv_mean and sv_std for the generic rand mat case
function svs = nonsymrandsv(N, M, xmean, xstd, Niter)
% get a set of singular values
Nsv = min([N,M]); % Number of singular values
svs = zeros(Niter,Nsv);
for iter=1:Niter
    H = xstd.*randn(N,M) + xmean;
    svs(iter,:) = svd(H);
end
