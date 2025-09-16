function MI = icatb_compute_mi(temp, temp2)
%% Compute FNC metrics using mutual information
%

temp  = icatb_remove_mean(temp');
temp = bsxfun(@rdivide, temp, std(temp));
%temp = temp*diag(1./std(temp));

if (exist('temp2', 'var'))
    temp2  = icatb_remove_mean(temp2');
    %temp2 = temp2*diag(1./std(temp2));
    temp2 = bsxfun(@rdivide, temp2, std(temp2));
end

order = size(temp, 2);

method = 1;
if (~exist('estpab', 'file'))
    method = 2;
end

if (~exist('temp2', 'var'))
    MI = zeros(order, order);
    for m = 1:order
        for n = m+1:order
            %MI(m, n) = sqrt(1-exp(-2*icatb_mutualinfo(temp(:, m), temp(:, n), method)));
            MI(m, n) = icatb_mutualinfo(temp(:, m), temp(:, n), method);
            MI(n, m) = MI(m, n);
        end
    end
else
    MI = zeros(size(temp, 2), size(temp2, 2));
    for m = 1:size(temp, 2)
        for n = 1:size(temp2, 2)
            %MI(m, n) = sqrt(1-exp(-2*icatb_mutualinfo(temp(:, m), temp2(:, n), method)));
            MI(m, n) = icatb_mutualinfo(temp(:, m), temp2(:, n), method);
        end
    end
end