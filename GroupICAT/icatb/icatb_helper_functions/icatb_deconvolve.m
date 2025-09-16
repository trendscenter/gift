function [y_recovered]= icatb_deconvolve(z, x)
% Deconvolution function

A = icatb_convmtx(x, length(z) - length(x) + 1);

y_recovered = pinv(A)*z;
