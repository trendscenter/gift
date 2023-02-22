% Description:
%
%     Create a design matrix for group-group interactions.
%
% Syntax:
%
%     [ X, terms ] = mGG2X(groups, offset)
%
% Inputs:
%
%     groups  - [ N x G ] (int)  - columns of qualitative variables
%     offset  - [ 1 x 1 ] (int)  - group offset index for terms
%     options - [ 1 x P ] (cell) - e.g. { 'verbose' }
%     columns - [ T x 2 ] (int)  - indices into the columns of groups
%
% Outputs:
%
%     X     - [ N x M ] (double)
%     terms - [ 1 x M ] (cell)
%
% Details:
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

function [ X, terms ] = mGG2X(groups, offset, options, columns)
    
    if ~exist('offset', 'var') || isempty(offset)
        offset = 0;
    end
    
    if ~exist('options', 'var')
        options = cell(0);
    end
    
    if ~exist('columns', 'var')
        columns = mNC2(size(groups, 2));
    end
    
    terms = {};
    X     = [];

    if ~isempty(strmatch('verbose', options, 'exact'))
        fprintf('\n')
    end
    
    for i = 1 : size(columns, 1)

        A = mG2X(groups(:, columns(i, 1)), offset);
        B = mG2X(groups(:, columns(i, 2)), offset);
        
        combinations = [];

        for j = 1 : size(A, 2)
            for k = 1 : size(B, 2)
                combinations(end + 1, :) = [ j k ];
            end
        end

        for j = 1 : size(combinations, 1)
            X = [ X A(:, combinations(j, 1)) .* B(:, combinations(j, 2)) ];
            terms{end + 1} = columns(i, :) + offset;
        end
        
        if ~isempty(strmatch('verbose', options, 'exact'))
            fprintf('Factor %d%d represents the interaction between columns %d and %d of groups and has %d levels.\n', ...
                columns(i, 1) + offset, columns(i, 2) + offset, ...
                columns(i, 1), columns(i, 2), size(combinations, 1))
        end
        
    end
    
