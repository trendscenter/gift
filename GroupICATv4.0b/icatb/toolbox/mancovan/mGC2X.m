% Description:
%
%     Create a design matrix for group-covariate interactions.
%
% Syntax:
%
%     [ X, terms ] = mGC2X(groups, covariates, a, b)
%
% Inputs:
%
%     groups     - [ N x G ] (int)    - columns of qualitative variables
%     covariates - [ N x C ] (double) - columns of quantitative variables
%     a          - [ 1 x 1 ] (int)    - group offset index for terms
%     b          - [ 1 x 1 ] (int)    - covariate offset index for terms
%     options    - [ 1 x P ] (cell)   - e.g. { 'verbose' }
%     columns    - [ T x 2 ] (int)    - see Details
%
% Outputs:
%
%     X     - [ N x M ] (double)
%     terms - [ 1 x M ] (cell)
%
% Details:
%
%     The first column of the columns input in this case is an index into the
%     columns of groups.  The second column of the columns input is an index
%     into the columns of covariates.  For example, [ 2 1 ] would refer to an
%     interaction between column 2 of groups and column 1 of covariates.
%
% Examples:
%
% Notes:
%
% Author(s):
%
%     William Gruner (williamgruner@gmail.com)
%
% References:
%
% Acknowledgements:
%
%     Many thanks to Dr. Erik Erhardt and Dr. Elena Allen of the Mind Research
%     Network (www.mrn.org) for their continued collaboration.
%
% Version:
%
%     $Author: williamgruner $
%     $Date: 2010-04-09 08:53:39 -0600 (Fri, 09 Apr 2010) $
%     $Revision: 491 $

function [ X, terms ] = mGC2X(groups, covariates, a, b, options, columns)
    
    if ~exist('a', 'var') || isempty(a)
        a = 0;
    end
    
    if ~exist('b', 'var') || isempty(b)
        b = 0;
    end
    
    if ~exist('options', 'var')
        options = {};
    end
    
    if ~exist('columns', 'var')
        
        columns = [];

        for i = 1 : size(groups, 2)
            for j = 1 : size(covariates, 2)
                columns(end + 1, :) = [ i j ];
            end
        end

    end
    
    terms = {};
    X     = [];
    
    [ G, g ] = mG2X(groups, 0);

    if ~isempty(strmatch('verbose', options, 'exact'))
        fprintf('\n')
    end
    
    for i = 1 : size(columns, 1)
        
        if ~isempty(strmatch('verbose', options, 'exact'))
            fprintf('Factor %d%d represents the interaction between column %d of groups and column %d of covariates.\n', ...
                columns(i, 1) + a, columns(i, 2) + b, columns(i, 1), columns(i, 2))
        end

        I = mFindTerms(columns(i, 1), g);
        x = G(:, I) .* repmat(covariates(:, columns(i, 2)), 1, length(I));
        X = cat(2, X, x);

        for j = 1 : length(I)
            terms{end + 1} = columns(i, :) + [ a b ];
        end
        
    end
