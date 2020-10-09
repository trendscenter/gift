function corr_coeff = icatb_corr2(x, y)
% computes correlation coefficient

meanX = mean(x(:));
meanY = mean(y(:));

% Remove mean
x = x - meanX;
y = y - meanY;

corr_coeff = sum(sum(x.*y))/sqrt(sum(sum(x.*x))*sum(sum(y.*y)));
