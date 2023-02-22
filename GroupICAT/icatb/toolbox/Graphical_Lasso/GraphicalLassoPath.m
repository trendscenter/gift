function [wList, thetaList, lambdaList, errors] = GraphicalLassoPath(pop, lambdaList, approximate, verbose, penalDiag, tolThreshold, maxIter)
% [wList, thetaList, lambdaList, errors] = GraphicalLassoPath(pop,
% lambdaList, approximate, verbose, penalDiag, tolThreshold, maxIter)
%
% Computes a list of regularized estimates of covariance matrix and its
% inverse along the regularization path
% Inputs:
%  - pop: the set of samples to be used for covariance estimation in an NxP
%  matrix where N is the number of samples and P is the number of variables
%  - lambdaList(o*): a vector containing different lambda values to be used
%  along the regularization path. If not specified a list of 10 increasing
%  lambda values will be generated using the entries of empirical
%  covariance matrix
%  - approximate(o): a flag indicating whether to use approximate
%  estimation (Meinhausen-Buhlmann approximation)
%  - verbose(o): a flag indicating whether to output algorithm process
%  - penalDiag(o): a flag indicating whether to penalize diagonal elements
%  of the covariance matrix
%  - tolThreshold(o): the amount of tolerance acceptable for covariance
%  matrix elements before terminating the algorithm
%  - maxIter(o): maximum number of iteration to perform in the algorithm
% *: o indicates optional arguments
% Outputs:
%  - wList: the list of estimations obtained for covariance matrix along
%  the regularization path
%  - thetaList: the list of estimations obtained for inverse covariance
%  matrix along the regularization path
%  - lambdaList: actual lambda values used along the regularization path
%  - errors: a list of error flags indicating the validity of each solution
%  along the regularization path
%
% Code by: Hossein Karshenas (hkarshenas@fi.upm.es)
% Date: 10 Feb 2011

if nargin < 1
    error('Too few input parameters.');
end
numVars = size(pop, 2);
sigma = cov(pop);
if nargin < 2
    tmp = max(max(abs(sigma)));
    lambdaList = (tmp ./ 10):(tmp ./ 10):tmp;
end
pathLength = length(lambdaList);
if nargin < 3
    approximate = 0;
end
if nargin < 4
    verbose = 0;
end
if nargin < 5
    penalDiag = 1;
end
if nargin < 6
    tolThreshold = 1e-4;
end
if nargin < 7
    maxIter = 1e4;
end

[wList, thetaList, lambdaList, errors] = glasso(numVars, sigma, 1, pathLength, lambdaList, approximate, verbose, penalDiag, tolThreshold, maxIter);
