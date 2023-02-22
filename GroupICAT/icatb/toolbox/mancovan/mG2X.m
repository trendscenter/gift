% Description:
%
%     Create a design matrix for groups.
%
% Syntax:
%
%     [ X, terms ] = mG2X(groups, offset)
%
% Inputs:
%
%     groups  - [ N x G ] (int)  - columns of qualitative variables
%     offset  - [ 1 x 1 ] (int)  - group offset index for terms
%     options - [ 1 x P ] (cell) - e.g. { 'verbose' }
%     columns - [ 1 x T ] (int)  - indices into columns of groups
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

function [ X, terms ] = mG2X(groups, offset, options, columns)
    
    if ~exist('offset', 'var') || isempty(offset)
        offset = 0;
    end
    
    if ~exist('options', 'var')
        options = cell(0);
    end
    
    if ~exist('columns', 'var')
        columns = 1 : size(groups, 2);
    end
    
    terms = {};
    X     = [];

    if ~isempty(strmatch('verbose', options, 'exact'))
        fprintf('\n')
    end

    for i = columns
        
        [ B, I, J ] = unique(groups(:, i));
        
        if ~isempty(strmatch('verbose', options, 'exact'))
            fprintf('Factor %d represents column %d of groups and has %d levels.\n', ...
                i + offset, i, length(B))
        end
        
%         for j = 1 : length(B) - 1
%             X(J ~=         j, length(terms) + j) =  0;
%             X(J ==         j, length(terms) + j) =  1;
%             X(J == length(B), length(terms) + j) = -1;
%         end
         
        for j = 2 : length(B) 
            X(J ~=         j, length(terms) + (j-1)) =  0;
            X(J ==         j, length(terms) + (j-1)) =  1;
        end



        for j = 1 : length(B) - 1
            terms{end + 1} = i + offset;
        end            
        
    end
    
