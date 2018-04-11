function mat = icatb_vec2mat(vec, full, fillNaN)
% mat = vec2mat(vec, full)
%returns the correlation matrix from a vector.
% If vec is [p x n*(n-1)/2], returns [p x n x n]
% if full == 1, returns the symmetric matrix

if nargin<2
    full = 1;
end

if (~exist('fillNaN', 'var'))
    fillNaN = 0;
end

if isvector(vec)
    N = length(vec);
    n = 1/2 + sqrt(1+8*N)/2;
    %mat = nan*ones(n);
    if (fillNaN)
        mat =  nan*ones(n);
    else
        mat = zeros(n);
    end
    temp = ones(n);
    %% find the indices of the lower triangle of the matrix
    IND = find((temp-triu(temp))>0);
    mat(IND) = vec;
    
    if full
        tempmat = flipud(rot90(mat));
        tempmat(IND) = vec;
        mat = flipud(rot90(tempmat));
    end
    
elseif ndims(vec) == 2
    [p,N] = size(vec);
    n = 1/2 + sqrt(1+8*N)/2;
    mat = zeros([p,n,n]);
    
    temp = ones(n);
    %% find the indices of the lower triangle of the matrix
    IND = find((temp-triu(temp))>0);
    
    for ii = 1:p
        tempmat  = zeros(n);
        tempmat(IND) = vec(ii,:);
        if full
            tempmat2 = flipud(rot90(tempmat));
            tempmat2(IND) = vec(ii,:);
            mat(ii,:,:) = flipud(rot90(tempmat2));
        else
            mat(ii,:,:) = tempmat;
        end
    end
    
end