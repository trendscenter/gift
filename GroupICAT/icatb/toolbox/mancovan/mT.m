% Description:
%
%     Compute t-statistics and p-values associated with pair-wise differences
%     between levels of a group or slopes of a regression line.
%
% Syntax:
%
%     [ t, p, stats ] = mT(Y, groups, [], term, options)
%     [ t, p, stats ] = mT(Y, [], covariates, term, options)
%     [ t, p, stats ] = mT(Y, groups, covariates, term, options)
%
%     [ t, p, stats ] = mT(Y, X, terms, term, options)
%
% Inputs:
%
%     Y          - [ N x M ] (double) - multivariate response
%     groups     - [ N x G ] (int)    - qualitative variables
%     covariates - [ N x C ] (double) - quantitative variables
%     term       - [ 1 x T ] (int)    - see Details
%     options    - [ 1 x P ] (cell)   - see Options
%
%     X     - [ N x T ] (double) - design matrix
%     terms - [ 1 x T ] (cell)   - model terms
%
% Outputs:
%
%     t - [ 1 x M ] (double)
%     p - [ 1 x M ] (double)
%
%     stats.Terms - [ 1 x B ] (cell) - full model terms numbering groups from
%         one through size(groups, 2) and covariates from size(groups, 2) + 1
%         through size(groups, 2) + size(covariates, 2) + 1, and interactions
%         according to combinations of these numbers.
%
%     stats.X - [ N x B ] (double) - full model design matrix.
%
%     stats.B - [ B x M ] (int) - full model regression coefficients associated
%         with each column of Y in the order presented in stats.Terms.
%
%     stats.SSE - [ 1 x M ] (int) - full model sum of squared errors associated
%         with each column of Y.
%
%     stats.DFE - [ 1 x M ] (int) - degrees of freedom used in the computation
%         of stats.MSE.
%
%     stats.MSE - [ 1 x M ] (int) - full model mean squared errors for each
%
%         column of Y.
%
%     stats.Term - [ 1 x T ] (int) - copy of the input term.
%
%     stats.Levels - [ L x T ] (int) - indices into stats.Terms indicating which
%         pair-wise differences were used for each t-test.  Rows of this output
%         correspond to rows of t and p and a full list of pair-wise differences
%         can be displayed by evaluating stats.Terms(stats.Levels).
%
% Details:
%
%     The term input may be a scalar index for a main effect, in which case this
%     input is an index into the columns of [ groups covariates ], or a 1 x 2
%     vector of indices, in which case this input is still an index into the
%     column-wise concatenation of groups and covariates, but with term(1)
%     indicating the first term in the interaction and term(2) indicating the
%     second term in the interaction.  For example, in a model with 2 groups and
%     2 covariates, the group-group interaction is given by term = [ 1 2 ], the
%     group-covariate interactions are given by term = [ 1 3 ], [ 1 4 ], [ 2 3 ],
%     and [ 2 4 ], and the covariate-covariate interaction is given by [ 3 4 ].
%
% Options:
%
%     'group-group'         - include group-group interactions.
%     'covariate-covariate' - include covariate-covariate interactions.
%     'group-covariate'     - include group-covariate interactions.
%     'verbose'             - display extra information to the command window.
%
% Examples:
%
%     The following example uses a simple additive model with covariates, but no
%     interactions and avoids using the Statistics Toolbox:
%
%         n          = 100;
%         groups     = round(3 * rand(n, 2) + 0.5);
%         covariates = 10 * randn(n, 2);
%         Y          = groups + covariates + randn(n, 2);
%
%     When the covariates are not included, the differences between levels of
%     the first group are insignificant:
%
%         [ t, p, stats ] = mT(Y, groups, [], 1)
%
%     When the covariates are included, the differences between levels of the
%     first group become apparent:
%
%         [ t, p, stats ] = mT(Y, groups, covariates, 1)
%
%     as do the differences between levels of the second group:
%
%         [ t, p, stats ] = mT(Y, groups, covariates, 2)
%
%     In addition, the significance of the slope associated with the first
%     covariate and each column of the response may be computed as follows:
%
%         [ t, p, stats ] = mT(Y, groups, covariates, 3)
%
%     and the slope associated with the second covariate and each column of the
%     response may be computed as follows:
%
%         [ t, p, stats ] = mT(Y, groups, covariates, 4)
%
%     An interaction may be introduced as follows:
%
%         Y = groups + covariates + [ groups(:, 1) .* covariates(:, 2) ...
%             zeros(n, 1) ] + randn(n, 2);
%
%     and evaluated as follows:
%
%         [ t, p, stats ] = mT(Y, groups, covariates, [ 1 4 ], ...
%             { 'group-covariate' 'verbose' })
%
%     As a second example, consider the more sophisticated multivariate response
%     given by:
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
%     The t-statistics and p-values associated with differences among levels of
%     the first group may be computed for each column of Y as follows:
%
%         [ t, p, stats ] = mT(Y, groups, covariates, 1, { 'verbose' });
%
%     The t-statistics and p-values associated with the slope of the regression
%     relationship between the first covariate and each column of Y may be
%     computed as follows:
%
%         [ t, p, stats ] = mT(Y, groups, covariates, 3 + 1, { 'verbose' });
%
%     and similarly for the remaining groups and covariates.  In each case, the
%     following plot reveals the location of the responsible component of L:
%
%         plot(-log10(p)')
%
% Notes:
%
% Author(s):
%
%     William Gruner (williamgruner@gmail.com)
%
% References:
%
%     Refer to the references listed in mancovan.m.
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

function [ t, p, stats ] = mT(Y, groups, covariates, term, options)

if nargin == 0
    t = BIT(); return
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

n = size(X, 1);
r = rank(X);
I = mFindTerms(term, terms);

if isempty(I)
    error('The specified term was not found in the model.')
end

% Create a matrix of indices for contrasts.

contrasts = [ I(:) , zeros(length(I), 1) ];

% Create a matrix of indices for pairs of contrasts.

combinations = mNC2(length(I));

% Concatenate all contrasts into a single matrix of indices.

contrasts = cat(1, contrasts, I(combinations));

% The full inverse is many times faster than the pseudoinverse.

if r == size(X, 2)
    invXTX = inv(X' * X);
else
    invXTX = pinv(X' * X);
end

M = X * invXTX * X';
invXTXXT = invXTX * X';

if ~isempty(strmatch('verbose', options, 'exact'))
    fprintf('\n')
    fprintf('Computing SSE ... ')
    tSSE = tic;
end

v = sscanf(version, '%d.%d.%d');
v = 10.^(0 : -1 : -(length(v) - 1)) * v;

if v > 7.5
    eyeM  = bsxfun(@minus, eye(size(M)), M);
    tempM = Y' * eyeM;
    SSE   = sum(bsxfun(@times, tempM', Y));
else
    eyeM  = eye(size(M)) - M;
    tempM = Y' * eyeM;
    SSE   = sum(tempM'.* Y);
end

if ~isempty(strmatch('verbose', options, 'exact'))
    fprintf('finished in %f seconds.\n', toc(tSSE))
end

B   = invXTXXT * Y;

DFE = n - r;
MSE = SSE / DFE;

if ~isempty(strmatch('verbose', options, 'exact'))
    fprintf('\n')
    fprintf('Computing t-statistics ...\n')
    fprintf('\n')
    tMT = tic;
end

t = repmat(NaN, size(contrasts,1), size(Y, 2));
p = repmat(NaN, size(contrasts,1), size(Y, 2));

for j = 1 : size(Y, 2)
    
    s = MSE(j) * invXTX;
    
    for i = 1 : size(contrasts, 1)
        
        if contrasts(i, 2) == 0
            
            % Test if a regression coefficient is different from zero.
            
            if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
                fprintf('\tT-Test %d for Y(:, %d): B(%d, %d) == 0 ... ', ...
                    i, j, contrasts(i, 1), j)
            end
            
            t(i, j) = B(contrasts(i, 1), j) / sqrt(s(contrasts(i, 1), contrasts(i,1)));
            p(i, j) = 2 * (1 - mTCDF(abs(t(i, j)), n - r));
            
            if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
                fprintf('p = %g\n', p(i, j));
            end
            
        else
            
            if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
                fprintf('\tT-Test %d for Y(:, %d): B(%d, %d) == B(%d, %d) ... ', ...
                    i, j, contrasts(i, 1), j, contrasts(i, 2), j)
            end
            
            % Test if two regression coefficinets are different from each other.
            
            t(i, j) = (B(contrasts(i, 1), j) - B(contrasts(i, 2), j)) / sqrt(s(contrasts(i, 1), contrasts(i, 1)) + ...
                s(contrasts(i, 2), contrasts(i, 2)) - 2 * s(contrasts(i, 1), contrasts(i, 2)));
            
            p(i, j) = 2 * (1 - mTCDF(abs(t(i, j)), n - r));
            
            if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
                fprintf('p = %g\n', p(i, j));
            end
            
        end
    end
    
    if ~isempty(strmatch('verbose', options, 'exact')) && j == 1
        fprintf('\n')
    end
    
end

levels = contrasts;
levels(levels > 0) = levels(levels > 0) - levels(1,1) + 1;

stats.Terms  = terms;
stats.X      = X;
stats.B      = B;
stats.SSE    = SSE;
stats.DFE    = DFE;
stats.MSE    = MSE;
stats.Term   = term;
stats.Levels = levels;

if ~isempty(strmatch('verbose', options, 'exact'))
    fprintf('... finished in %f seconds.\n', toc(tMT))
    fprintf('\n')
end

function b = BIT()

b = true;

% Compare the results to those obtained from the previous version.

rand ('seed', 0)
randn('seed', 0)

n          = 1000;
groups     = round(3 * rand(n, 2) + 0.5);
covariates = randn(n, 2);

Y = groups + covariates + ...
    [ groups(:, 1) .* groups(:, 2) zeros(n, 1) ] + randn(n, 2);

[ t, p, stats ] = mT(Y, groups, covariates, 1, ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [
    -63.659577493314664  -1.558037833198396
    -1.429399383538909   0.323432588359129
    -36.597375547607541  -1.099416539613320 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [
    0   0.119547139350506
    0.153208160657149   0.746436656117869
    0   0.271856423543536 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, 2, ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [
    -43.998644412972560 -23.374016759590667
    0.037663118597656  -0.345362811487873
    -25.331536619358733 -13.247269701786575 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [
    0                   0
    0.969963952692249   0.729895746141207
    0                   0 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, 3, ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [ 32.426050257963418 0.814869311153052 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [ 0 0.415345098645853 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, 4, ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [ 0.190124598399025 33.745828629822775 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [ 0.849250868718951 0 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, [ 1 2 ], ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [
    16.825108617220536   1.213285095215801
    -1.072528078927162  -0.745163499279065
    -1.484147052326498   0.383450655931146
    0.348597922344957   1.715706947088990
    10.230585919645254   1.119198885949145
    10.886542515174956   0.508970824908132
    13.831421607018257  -0.371205671521071
    0.278671923770819  -0.943270442443392
    -0.821539781478530  -1.400465391034129
    -1.051239470349919  -0.784428263560105 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [
    0   0.225313189236329
    0.283746944057963   0.456351554082997
    0.138091460795758   0.701468816848293
    0.727466062780892   0.086531760982572
    0   0.263329539424501
    0   0.610887243735730
    0   0.710564581350285
    0.780555401217813   0.345775068248919
    0.411538687453927   0.161690527387597
    0.293407765441479   0.432978425049231 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, [ 3 4 ], ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [ -0.031159942230402 -0.810374331050160 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [ 0.975148330481621 0.417921804493396 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, [ 1 3 ], ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [
    -0.172948518741838  -0.122183498481013
    -1.401051585421879   0.807441149728400
    0.690725036516936  -0.529168338935551 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [
    0.862727606180638   0.902778744641907
    0.161515223042944   0.419608296451506
    0.489901915528637   0.596808514224463 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, [ 1 4 ], ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [
    -0.335066004029774   1.269700725980901
    1.690743420205252  -0.359649120709498
    -1.171203576567804   0.939716170495065 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [
    0.737647002754556   0.204492667192531
    0.091203798888301   0.719187029037324
    0.241801643694997   0.347594723721451 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, [ 2 3 ], ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [
    1.172807004587313   0.528898804990101
    -0.508202969391037   0.127967553999966
    0.967439760395387   0.232598373964932 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [
    0.241158199656773   0.596995416675322
    0.611425378979106   0.898200904085521
    0.333562926466858   0.816121856651645 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

[ t, p, stats ] = mT(Y, groups, covariates, [ 2 4 ], ...
    { 'group-group' 'covariate-covariate' 'group-covariate' 'verbose' });

e = t - [
    0.361860128953282   0.036967033566348
    -1.179649790281506  -0.112984082460597
    0.892341502739559   0.086797959331714 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end

e = p - [
    0.717534564266703   0.970518819684806
    0.238425814127398   0.910066314162987
    0.372428967150598   0.930849840649059 ];

if any(abs(e(:)) > sqrt(eps))
    b = false;
end


