function  y = icatb_lap_pdf(x,a, b)
% y = icatb_Lap_pdf(x,a, b)
%
% Created By: Baoming Hong, IOL, CT, USA
% Data: 1-08-04
% Laplace distribution f(x)= 1/(2b)*exp(-|x-a|/b). N is the random data number 
% you want to create.  "b" is a shape parameter (b>0); "a" is a shift factor.
% the expected value: E(X) = a; its variance: Var(X)= 2*b^2;, skewness=0,
% kurtosis = 3;

 z = (abs(x-a))/b;
 y = exp(-z)/(2*b);
 
 % --------- Stephen Hong, 4/25/03 -----------------
% function  Lappdf = Lap_pdf(a, b, N) 
%  i = 0;
%  for X=-N:0.1:N;
%     i= i + 1;
%     Y= (abs(X-a))/b;
%     Lappdf(i) = exp(-Y)/(2*b);  % 
%  end