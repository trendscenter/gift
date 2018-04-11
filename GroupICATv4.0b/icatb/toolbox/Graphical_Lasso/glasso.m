function [w, theta, iter, avgTol, hasError] = glasso(numVars, s, computePath, lambda, approximate, warmInit, verbose, penalDiag, tolThreshold, maxIter, w, theta)
% [w, theta, iter, avgTol, hasError] = glasso(numVars, s, computePath,
% lambda, approximate, warmInit, verbose, penalDiag, tolThreshold, maxIter,
% w, theta)
%
% OR
%
% [wList, thetaList, lambdaList, errors] = glasso(numVars, s, computePath,
% pathLength, lambdaList, approximate, verbose, penalDiag, tolThreshold,
% maxIter)
%
% The Matlab's mex interface function for Graphical Lasso algorithm which
% obtains a regularized estimate of covariance matrix and its inverse
% It has two versions: 1) obtaining a single solution for a given lambda
%                      2) obtaining a list of solutions along the
%                      regularization path
% Input:
%  - numVars: number of problem variables
%  - s: the empirical computation of covariance matrix
%   computePath: a flag indicating whether to compute the regularization
%   path or not. This flag specifies how the following parameters are going
%   to be treated.
%  - lambda: the matrix of regularization penalties
%  - approximate: a flag indicating whether to use approximate estimation
%  (Meinhausen-Buhlmann approximation)
%  - warmInit: a flag indicating whether the estimation will start from
%  given initial values of 'w' and 'theta'
%  - verbose: a flag indicating whether to output algorithm process
%  - penalDiag: a flag indicating whether to penalize diagonal elements of
%  the covariance matrix
%  - tolThreshold: the amount of tolerance acceptable for covariance matrix
%  elements before terminating the algorithm
%  - maxIter: maximum number of iteration to perform in the algorithm
%  - w: the initial value of covariance matrix used for warm initialization
%  - theta: the initial value of inverse covariance matrix used for warm
%  initialization
%  - pathLength: the number of different lambda values along the
%  regularization path
%  - lambdaList: different lambda values along the regularization path
% Output:
%  - w: the estimated covariance matrix
%  - theta: the estimated inverse covariance matrix
%  - iter: actual number of iterations performed in the algorithm
%  - avgTol: average tolerance of covariance matrix entries before
%  terminating the algorithm
%  - hasError: a flag indicating whether the algorithm terminated
%  erroneously or not
%  - wList: the list of estimations obtained for covariance matrix along
%  the regularization path
%  - thetaList: the list of estimations obtained for inverse covariance
%  matrix along the regularization path
%  - lambdaList: actual lambda values used along the regularization path
%  - errors: a list of error flags indicating the validity of each solution
%  along the regularization path
%
% The code is adapted for Matlab by Hossein Karshenas (hkarshenas@fi.upm.es)
% from the original code provided by authors of the following paper:
% Jerome Friedman, Trevor Hastie and Robert Tibshirani (2007).
% Sparse inverse covariance estimation with the lasso.
% Biostatistics 2007. http://www-stat.stanford.edu/~tibs/glasso/
%
% It also uses some of the utility functions provided by Hui Jiang
% (jiangh@stanford.edu) in his Matlab wrapper written for generalized
% linear models estimation: http://www-stat.stanford.edu/~tibs/glmnet-matlab/
%
% License: GPL, 10 Fubruary 2011

error('mex file absent, type ''mex glasso.F'' to compile');
