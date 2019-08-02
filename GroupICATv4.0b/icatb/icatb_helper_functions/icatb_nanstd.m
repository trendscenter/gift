function  s = icatb_nanstd(x, dim)
%% Compute nanstd
%

if (numel(x) == length(x))
    x = x(:);
end


if (~exist('dim', 'var'))
    dim = 1;
end

sz = size(x);
if dim > length(sz)
    sz(end+1:dim) = 1;
end
fac = ones(size(sz));
fac(dim) = sz(dim);

nans = isnan(x);

%% Center data
m = icatb_nanmean(x, dim);
x = (x - repmat(m, fac)).^2;

%% Std
n = sum(~nans, dim)+eps;
n = max(1, n - 1);
x(nans) = 0;
x = x.*(~nans);
s = sqrt(sum(x, dim)./n);