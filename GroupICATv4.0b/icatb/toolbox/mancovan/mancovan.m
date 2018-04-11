% Description:
%
%     Compute the Lawley-Hotelling trace and the associated p-value for each
%     term in a MANCOVAN model.
%
% Syntax:
%
%     [ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, groups)
%     [ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, [], covariates)
%     [ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, groups, covariates)
%     [ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, groups, covariates, options)
%
%     [ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, X, terms, options)
%
% Inputs:
%
%     Y          - [ N x M ] (double) - multivariate response
%     groups     - [ N x G ] (int)    - qualitative variables
%     covariates - [ N x C ] (double) - quantitative variables
%     options    - [ 1 x P ] (cell)   - see Options
%
%     X     - [ N x T ] (double) - design matrix
%     terms - [ 1 x T ] (cell)   - model terms
%
% Outputs:
%
%     T        - [ (G + C) x 1 ] (double)
%     p        - [ (G + C) x 1 ] (double)
%     FANCOVAN - [ (G + C) x M ] (double)
%     pANCOVAN - [ (G + C) x M ] (double)
%
%     stats.U - [ N x N ] (double) - U resulting from the SVD (option 'SVD').
%     stats.S - [ N x N ] (double) - S resulting from the SVD (option 'SVD').
%     stats.V - [ M x N ] (double) - V resulting from the SVD (option 'SVD').
%
%     stats.BIC - [ 1 x P ] (double) - values of the Bayesian information
%         criterion (BIC) associated with the SVD (option 'SVD').
%
%     stats.PVE - [ 1 x N ] (double) - the percentage of variance in Y explained
%         by each consecutive column of U (option 'SVD').
%
%     stats.Terms - [ 1 x B ] (cell) - full model terms numbering groups from
%         one through size(groups, 2) and covariates from size(groups, 2) + 1
%         through size(groups, 2) + size(covariates, 2) + 1, and interactions
%         according to combinations of these numbers.
%
%     stats.X - [ N x B ] (double) - full model design matrix, including columns
%         for the grand mean, groups, covariates, and interactions.
%
%     stats.Y - [ N x M ] (double) - the multivariate response used for the
%         computation.  This will be equal to U if an SVD was performed.
%
%     stats.B - [ B x M ] (double) - regression coefficients associated with the
%         full model.
%
%     stats.SSE - [ 1 x M ] (int) - full model sum of squared errors associated
%         with each column of Y.
%
%     stats.DFE - [ 1 x M ] (int) - degrees of freedom used in the computation
%         of stats.MSE.
%
%     stats.MSE - [ 1 x M ] (int) - full model mean squared errors associated
%         with each column of Y.
%
% Details:
%
%     In addition to the MANCOVAN results, F-statistics and the associated
%     p-values are computed for each group, each covariate, and (optionally)
%     each two-way interaction and for each column of Y using ANCOVAN models.
%     These can be useful in isolating sources of significance found in the
%     MANCOVAN model.
%
% Options:
%
%     'group-group'         - include group-group interactions.
%     'covariate-covariate' - include covariate-covariate interactions.
%     'group-covariate'     - include group-covariate interactions.
%     'SVD'                 - reduce the dimensionality of Y using an SVD.
%     'verbose'             - display extra information to the command window.
%
% Examples:
%
%     The following example uses a simple additive model with covariates, but no
%     interactions and avoids using the Statistics Toolbox:
%
%         n          = 100; 
%         groups     = round(2 * rand(n, 2) + 0.5);
%         covariates = 10 * randn(n, 2);
%         Y          = groups + covariates + randn(n, 2);
%
%         [ T, p ] = mancovan(Y, groups, covariates, { 'verbose' });
%
%     Depending on the randomization, this should result in highly-significant
%     group effects.  The same is not true when covariates are excluded:
%
%         [ T, p ] = mancovan(Y, groups, []);
% 
%     The following example makes use of the SVD option:
%
%         n = 100;
%         m = 6000;
%         x = -3 : 6 / (m - 1) : 3;
%         z = zeros(n, 5);
%         e = randn(n, 13);
%
%         groups     = round(3 * rand(n, 3) + 0.5);
%         covariates = randn(n, 3);
%
%         L = normpdf(repmat(x, 13, 1), repmat([ -3 : 0.5 : 3 ]', 1, m), 0.25);
%         Y = ([ z groups z ] + [ z covariates z ] + e) * L + randn(n, m) / 20;
%
%         [ T, p ] = mancovan(Y, groups, covariates, { 'SVD' 'verbose' })
%
%     To test the same example when the effects are purely random, replace the
%     assignment to Y above with:
%
%         Y = randn(n, 13) * L + 0.05 * randn(n, m);
%
%     To introduce a significant two-way interaction between the second and
%     third groups, replace the assignment to Y above with:
%
%         Y = ([ z groups(:, 1:3) groups(:, 2) .* groups(:, 3) z(:, 1:4) ] + ...
%             randn(n, 13)) * L + randn(n, m) / 20;
%
%     and use the 'group-group' option to test for two-way interactions:
%
%         [ T, p ] = mancovan(Y, groups, covariates, { 'SVD' 'group-group' })
%
% Notes:
%
% Author(s):
%
%     William Gruner (williamgruner@gmail.com)
%
% References:
%
%     [1] Advanced Linear Modeling, Second Edition, Ronald Christensen
%     [2] Applied Linear Statistical Models, Fourth Edition, Neter et. al.
%     [3] http://en.wikipedia.org/wiki/Bayesian_information_criterion
%     [4] http://en.wikipedia.org/wiki/F-distribution
%     [5] http://en.wikipedia.org/wiki/Student's_t-distribution
%     [6] http://www.statsoft.com/textbook/general-linear-models
%
% Acknowledgements:
%
%     Many thanks to Dr. Erik Erhardt and Dr. Elena Allen of the Mind Research
%     Network (www.mrn.org) for their continued collaboration.
%
% Version:
%
%     $Author: williamgruner $
%     $Date: 2010-04-09 10:23:30 -0600 (Fri, 09 Apr 2010) $
%     $Revision: 492 $

function [ T, p, FANCOVAN, pANCOVAN, stats ] = mancovan(Y, groups, covariates, options)
    
    if nargin == 0
        T = BIT(); return
    end
    
    if ~exist('covariates', 'var')
        covariates = [];
    end
    
    if ~exist('options', 'var')
        options = cell(0);
    end
    
    if iscell(covariates)
        X     = groups;
        terms = covariates;
    else
        [ X, terms ] = mX(groups, covariates, options);
    end
        
    if ~isempty(strmatch('SVD', options, 'exact'))
        
        [ BIC, U, S, V ] = mSVD(Y, options);

        stats.U   = U;
        stats.S   = S;
        stats.V   = V;
        stats.BIC = BIC;
        stats.PVE = cumsum(diag(S) ./ sum(diag(S)));
        
        b = find(diff(BIC) > 0);
        b = b(1);
        U = U(  :, 1:b);
        S = S(1:b, 1:b);
        V = V(  :, 1:b);
        
    else
        
        U = Y;
        
    end
    
    if size(U, 2) > size(U, 1) - rank(X) - 2
        error(sprintf('MANCOVAN needs %d more degree(s) of freedom to proceed.  You may increase the number of samples, decrease the number of terms, or decrease the dimensionality of Y using the SVD option.', ...
            size(U, 2) - (size(U, 1) - rank(X) - 2)))
    end
    
    B        = mUnique(terms(2:end));
    M        = X * pinv(X' * X) * X';
    T        = repmat(NaN, length(B), 1);
    p        = repmat(NaN, length(B), 1);
    FANCOVAN = repmat(NaN, length(B), size(U, 2));
    pANCOVAN = repmat(NaN, length(B), size(U, 2));

    if ~isempty(strmatch('verbose', options, 'exact'))
        fprintf('\n')
    end
            
    for i = 1 : length(B)
        
        I0 = mTerms(B{i}, terms);
        X0 = X(:, I0);
        M0 = X0 * pinv(X0' * X0) * X0';
        
        if ~isempty(strmatch('verbose', options, 'exact'))
            mDispModels([ {0} B ], mUnique(terms(I0)))
        end
        
        [ T(i), p(i) ] = mLHT(U, X, X0, M, M0, options);
        
        for j = 1 : size(U, 2)
            [ FANCOVAN(i, j), pANCOVAN(i, j) ] = mF(U(:, j), X, X0, M, M0);
        end
        
        if ~isempty(strmatch('verbose', options, 'exact'))
            fprintf('\n')
        end
        
    end
    
    stats.Terms = terms;
    stats.X     = X;
    stats.Y     = Y;
    stats.B     = pinv(X' * X) * X' * U;
    stats.SSE   = U' * (eye(size(M)) - M) * U;
    stats.DFE   = size(X, 1) - rank(X);
    stats.MSE   = stats.SSE / stats.DFE;
    
function b = BIT()
    
    b = true;
    
    % Compare the results to those obtained using ANOVAN.

    rand ('seed', 0)
    randn('seed', 0)
    
    n      = 20;
    groups = round(3 * rand(n, 2) + 0.5);
    Y      = groups + randn(n, 2);
    
    [ T, pT, F, pF, stats ] = mancovan(Y(:, 1), groups, [], { 'verbose' });
    
    [ P, T, STATS ] = anovan(Y(:, 1), { groups(:, 1) groups(:, 2) }, ...
        'display', 'off');
    
    % STATS.coeffs(4) = -(stats.B(2) + stats.B(3))
    % STATS.coeffs(7) = -(stats.B(4) + stats.B(5))
    
    e = [
        stats.B(1) - STATS.coeffs(1)
        stats.B(2) - STATS.coeffs(2)
        stats.B(3) - STATS.coeffs(3)
        stats.B(4) - STATS.coeffs(5)
        stats.B(5) - STATS.coeffs(6) ];
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end
    
    e = P - pF;
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end
    
    [ T, pT, F, pF, stats ] = mancovan(Y(:, 2), groups, [], { 'verbose' });
    
    [ P, T, STATS ] = anovan(Y(:, 2), { groups(:, 1) groups(:, 2) }, ...
        'display', 'off');
    
    % STATS.coeffs(4) = -(stats.B(2) + stats.B(3))
    % STATS.coeffs(7) = -(stats.B(4) + stats.B(5))
    
    e = [
        stats.B(1) - STATS.coeffs(1)
        stats.B(2) - STATS.coeffs(2)
        stats.B(3) - STATS.coeffs(3)
        stats.B(4) - STATS.coeffs(5)
        stats.B(5) - STATS.coeffs(6) ];
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end
    
    e = P - pF;
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end
    
    % Compare the results to those obtained from the previous version.
    
    [ T, pT, F, pF, stats ] = mancovan(Y, groups, [], { 'verbose' });
    
    e = stats.B - [
         2.183165187602118   2.035473809855300
        -0.391354062627052  -0.257117252131537
        -0.450424640193000  -0.202241625716015
        -0.038453801019392  -0.861629718426372
         0.128166317053669   0.162254725283209 ];
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end
    
    e = T - [
         7.589356026429847
        11.028135346312721 ];
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end

    e = pT - [
        0.191544857647959
        0.082196744321958 ];
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end

    e = F - [
        3.479859685613939   1.252230405151584
        0.075170700614001   4.790685342519653 ];
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end

    e = pF - [
        0.057343414455294   0.314100704058015
        0.927932310615466   0.024611536245152 ];
    
    if any(abs(e(:)) > sqrt(eps))
        b = false;
    end


  
    