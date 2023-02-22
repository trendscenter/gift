function  m = icatb_nanmean(x, dim)
%% Compute nanmean

if (~exist('dim', 'var'))
    dim = 1;
end

try
    
    m = mean(x, dim, 'omitnan');
    
catch
    
    if (numel(x) == length(x))
        x = x(:);
    end
    
    nans = isnan(x);
    x(nans) = 0;
    
    n = sum(~nans, dim)+eps;
    x = x.*(~nans);
    m = sum(x, dim) ./ n;
    
end