function [Bhat, R2] = icatb_ls_fit(X,y)
%function [Bhat, R2] = ls_fit(X,y)
%does least squares minimization
%X = vector of explanatory variables (model)
%y = data vector
%Bhat = LS estimates of parameters (coefficients)
%R2 correlation coefficient

if ((size(X,2)>size(X,1))&(size(X,2)==size(y,2))),
   X = X';
   y=y';
end;

%Bhat = inv(X'*X)*X'*y;
% solve system of equations instead of inverting the matrix
Bhat = (X'*X) \ (X'*y);

N =length(y);
m2 = mean(y)*mean(y);
R2 = abs((Bhat'*X'*y - N*m2)/(y'*y-N*m2));
