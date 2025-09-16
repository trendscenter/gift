function y = icatb_logistic_pdf(x,m,b)
% y = icatb_logistic_pdf(x,m,b)
%
% Created By: Baoming Hong, IOL, CT, USA
% Data: 1-08-04
z = (x-m)/b;
y =1/b*exp(-z)./(1 + exp(-z)).^2;

% E(x) = m;
% Var(x) =1/3*pi^2*b^2;
% skewness = 0;
% kurtosis = 1.2;
