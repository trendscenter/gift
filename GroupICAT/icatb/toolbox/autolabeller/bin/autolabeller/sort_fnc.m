% This function modularizes an FNC matrix
% Inputs: 
%     fnc:     unsorted FNC matrix
%     labels:  A 3-column cell. The columns are:
%         IC#:                     should correspond to the row/column number of the FNC matrix input
%         network/noise label:     a column of 1/0 in case the FNC matrix contains both network and noise entries
%         functional domain label: a string containing the functional domain name of the IC. ICs belonging to each domain is modularized separately.
%         
%        Example labels input:
%
%         35,1,Visual
%         16,1,Default
%         36,1,Subcortical
%         ...
% 
% Prerequisite:
%     The reorder_mod function comes from the BCT toolbox- https://github.com/brainlife/BCT/tree/main/BCT/2019_03_03_BCT. It needs to be in Matlab path.
% 
% Note:
%     If the whole matrix needs to be modularized and no label is needed, then a dummy labels cell like this should work-
%         
%         1,1,""
%         2,1,""
%         3,1,""
%         ...
% 
function [network_idx_reordered, reordered_matrix, order_] = sort_fnc( fnc, labels )
    % sorted domain labels
    network_idx = find( cell2mat( labels(:,2) ) );
    noise_idx = find( ~cell2mat( labels(:,2) ) );
    network_fnc = fnc( network_idx, network_idx );

    % reorder modules
    [order_, reordered_matrix] = reorder_mod( network_fnc, labels( network_idx, 3 ) );

    % sorted output index of the networks
    network_idx_reordered = network_idx( order_ );

    % sorted output index of all components
    order_ = [network_idx_reordered; noise_idx];

    disp('done reordering FNC')

