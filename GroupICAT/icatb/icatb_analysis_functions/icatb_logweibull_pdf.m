function  y = icatb_logweibull_pdf(x,a,b)
% y = icatb_logweibull_pdf(x,a,b)
% Created By: Baoming Hong, IOL, CT, USA
% Data: 1-08-04

z = (a-x)/b;
y = 1/b*exp(z).*exp(-exp(z));
% parameter: a is shift factor, b is shape factor;
% E(X)= a+b*gamma;   Eluer-mascheroni constant: 0.5772156649
% Var(X)=1/6*pi^2*b^2;
% Kurtosis(X)= 12/5;