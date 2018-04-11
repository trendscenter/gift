function MI = icatb_compute_mi(temp)
%% Compute FNC metrics using mutual information
%

temp  = icatb_remove_mean(temp');
temp = temp*diag(1./std(temp));

order = size(temp, 2);
MI = zeros(order, order);

method = 1;
if (~exist('estpab', 'file'))
    method = 2;
end

for m = 1:order
    for n = m+1:order
        MI(m, n) = sqrt(1-exp(-2*icatb_mutualinfo(temp(:, m), temp(:, n), method)));
        MI(n, m) = MI(m, n);
    end
end