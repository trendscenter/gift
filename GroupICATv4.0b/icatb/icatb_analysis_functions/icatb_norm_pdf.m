function y = icatb_norm_pdf(x, mu, sigma)
%% Normal probability distribution function (pdf)
%
% Inputs:
% 1. x - Data
% 2. mu - Mean
% 3. sigma - Standard deviation
%
% Outputs:
% y - Normal pdf

if (~exist('mu', 'var'))
    mu = 0;
end

if (~exist('sigma', 'var'))
    sigma = 1;
end

y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);