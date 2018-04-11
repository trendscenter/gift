function [MN, low_inds] = icatb_domain_avg_fnc(M, compNStruct)
%% Average cells of FNC matrix
%
% Inputs:
% 1. M - Matrix of correlations
% 2. compNStruct - Data structure containing fields name and value.
%
% Outputs:
% MN - Network averaged correlations
%


for nF = 1:size(M, 1)
    tmp = compute_fnc(M(nF, :), compNStruct);
    if (nF == 1)        
        II = ones(size(tmp));
        II2 = tril(II);
        II2(isnan(tmp) == 1)=0;
        low_inds = find(II2 == 1);
        MN = zeros(size(M, 1), length(low_inds));
    end
    MN(nF, :) = tmp(low_inds);
end


function MN = compute_fnc(M, compNStruct)

if (numel(M) == length(M))
    M = icatb_vec2mat(M);
end

inds = (eye(size(M)))==1;
M(inds) = NaN;

MN = zeros(length(compNStruct), length(compNStruct));

endInds = 0;
for n = 1:length(compNStruct)
    
    startInds = endInds + 1;
    endInds = endInds + length(compNStruct(n).value);
    
    rowEndInds = 0;
    for m = 1:length(compNStruct)
        
        rowStartInds = rowEndInds + 1;
        rowEndInds = rowEndInds + length(compNStruct(m).value);
        subM = M(startInds:endInds, rowStartInds:rowEndInds);
        
        tmpMN = icatb_nanmean(subM(:));
        
        if (tmpMN == 0)
            tmpMN = NaN;
        end
        
        MN(n, m) = tmpMN;
        
    end
    
end

%MN = icatb_mat2vec(MN);
%MN = MN(:)';
