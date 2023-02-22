function [w, theta, iter, avgTol, hasError] = GraphicalLasso(pop, lambda, initStruct, approximate, warmInit, verbose, penalDiag, tolThreshold, maxIter, w, theta)
% [w, theta, iter, avgTol, hasError] = GraphicalLasso(pop, lambda,
% initStruct, approximate, warmInit, verbose, penalDiag, tolThreshold,
% maxIter, w, theta)
%
% Computes a regularized estimate of covariance matrix and its inverse
% Inputs:
%  - pop: the set of samples to be used for covariance estimation in an NxP
%  matrix where N is the number of samples and P is the number of variables
%  - lambda: the regularization penalty. Can be a single number or a matrix
%  of penalization values for each entry in the covariance matrix
%  - initStruct(o*): a matrix of size PxP, where zero entries will force
%  the corresponding entries in the inverse covariance matrix to be zero.
%  - approximate(o): a flag indicating whether to use approximate estimation
%  (Meinhausen-Buhlmann approximation)
%  - warmInit(o): a flag indicating whether the estimation will start from
%  given initial values of 'w' and 'theta'
%  - verbose(o): a flag indicating whether to output algorithm process
%  - penalDiag(o): a flag indicating whether to penalize diagonal elements of
%  the covariance matrix
%  - tolThreshold(o): the amount of tolerance acceptable for covariance matrix
%  elements before terminating the algorithm
%  - maxIter(o): maximum number of iteration to perform in the algorithm
%  - w(o): the initial value of covariance matrix used for warm initialization
%  - theta(o): the initial value of inverse covariance matrix used for warm
%  initialization
% *: o indicates optional arguments
% Outputs:
%  - w: the estimated covariance matrix
%  - theta: the estimated inverse covariance matrix
%  - iter: actual number of iterations performed in the algorithm
%  - avgTol: average tolerance of covariance matrix entries before
%  terminating the algorithm
%  - hasError: a flag indicating whether the algorithm terminated
%  erroneously or not
%
% Code by: Hossein Karshenas (hkarshenas@fi.upm.es)
% Date: 10 Feb 2011

if nargin < 2
    error('Too few input parameters.');
end
numVars = size(pop, 2);
[m, n] = size(lambda);
if m ~= n
    error('Regularization coefficients matrix should be symmetric matrix.');
elseif m > 1 && m ~= numVars
    error('Regularization coefficients matrix should have a size equal to the number of variables.');
end
if m == 1
    lambda = lambda .* ones(numVars);
end
if nargin > 2 && ~isempty(initStruct)
    initStruct = 1 - initStruct;
    initStruct = 10e9 .* initStruct;
    lambda = lambda + initStruct;
end
if nargin < 4
    approximate = 0;
end
if nargin < 5
    warmInit = 0;
end
if nargin < 6
    verbose = 0;
end
if nargin < 7
    penalDiag = 1;
end
if nargin < 8
    tolThreshold = 1e-4;
end
if nargin < 9
    maxIter = 1e4;
end
if nargin < 10
    if warmInit
        error('In warm initialization mode starting values for the covariance and precision matrices should be determined.');
    else
        w = zeros(numVars);
    end
end
if nargin < 11
    if warmInit
        error('In warm initialization mode starting values for the precision matrix should be determined.');
    else
        theta = zeros(numVars);
    end
end

[w, theta, iter, avgTol, hasError] = glasso(numVars, cov(pop), 0, lambda, approximate, warmInit, verbose, penalDiag, tolThreshold, maxIter, w, theta);

if hasError
    warning('The execution of the algorithm caused errors');
end
