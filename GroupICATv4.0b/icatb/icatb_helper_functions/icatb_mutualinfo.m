function M = icatb_mutualinfo(a, b, method)
%% Compute mutual info
%
% Inputs:
% a and b are variables of the same length

if (method == 1)
    M = mutualinfo(a, b);
else
    M = doMI(a, b);
end


function M = doMI(a, b)
%% Mutual information: use 2d histogram
%

bins1 = round(max(a) - min(a) + 1);
bins2 = round(max(b) - min(b) + 1);
xvert = linspace(min(a), max(a), bins1);
yvert = linspace(min(b), max(b), bins2);
h = icatb_histnd([a(:), b(:)], xvert, yvert);
h = h./sum(h(:)); % normalized joint histogram
py = sum(h, 1);
px = sum(h, 2);

% Mutual information
M = entropy(py) + entropy(px) - entropy(h(:));

function Hy = entropy(py)
%% Compute entropy
%

Hy = py.*log2(py);
Hy(isfinite(Hy)==0) = 0;
Hy = -sum(Hy);