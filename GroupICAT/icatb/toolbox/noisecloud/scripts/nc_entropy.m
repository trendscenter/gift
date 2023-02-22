function Hy = nc_entropy(a)
%% Compute entropy
%

a = double(a);
a = a(:);
bins1 = round(max(a) - min(a) + 1);
xvert = linspace(min(a), max(a), bins1);
py = hist(a, xvert);
py = py./sum(py(:)); 
Hy = double(py).*log2(double(py));
Hy(isfinite(Hy)==0) = 0;
Hy = -sum(Hy);