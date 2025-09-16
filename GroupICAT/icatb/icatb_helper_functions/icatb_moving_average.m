function Y = icatb_moving_average(X,F)
%MOVING_AVERAGE   Smooths a vector through the moving average method.
%   Y = MOVING_AVERAGE(X,F) Quickly smooths the vector X via averaging each 
%   element with the F elements at his right and the F elements at his 
%   left. The elements at the ends are also averaged but the extrems are
%   left intact.
%
%   Example:
%      x = 2*pi*linspace(-1,1); 
%      yn = cos(x) + 0.25 - 0.5*rand(size(x)); 
%      ys = nanmoving_average(yn,4,'interpnan');
%      plot(x,yn,x,ys), legend('noisy','smooth',4), axis tight
%
%   See also MOVING_AVERAGE2, NANMOVING_AVERAGE, NANMOVING_AVERAGE2

%   Written by
%   Lic. on Physics Carlos Adrián Vargas Aguilera
%   Physical Oceanography MS candidate
%   UNIVERSIDAD DE GUADALAJARA 
%   Mexico,  march 2006
%
%   nubeobscura@hotmail.com

% Moving average method, except the ends:
Y = boxcar_window(X,F);

% Smooths the ends:
Wwidth = 2*F+1;
N = length(X);
% First end:          
Yini = cumsum(X(1:Wwidth-2));     
Yini = Yini(:).';
Yini = Yini(1:2:end)./(1:2:Wwidth-2);            
Y(1:F) = Yini;
% Last end:
Yfin = cumsum(X(N:-1:N-Wwidth+3));
Yfin = Yfin(:).';
Yfin = Yfin(end:-2:1)./(Wwidth-2:-2:1);
Y(N-F+1:N) = Yfin;


function Y = boxcar_window(X,F)
% Boxcar window of length 2F+1 via recursive moving average (really fast)
%
% nubeobscura@hotmail.com


if F == 0
 Y = X;
 return
end

N = length(X);
Wwidth = 2*F + 1;                     % filter width
Y = zeros(size(X));          
Y(F+1) = sum(X(1:Wwidth));    
for n = F+2:N-F
 Y(n) = Y(n-1) + X(n+F) - X(n-F-1);   % recursive moving average
end
Y = Y/Wwidth;


% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com