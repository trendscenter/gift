function Y = icatb_nanmoving_average(X,F,varargin)
%NANMOVING_AVERAGE   Smooths a vector through the moving average method.
%   Y = NANMOVING_AVERAGE(X,F) Quickly smooths the vector X via averaging 
%   each element with the F elements at his right and the F elements at his 
%   left, ignoring NaN's. The elements at the ends are also averaged but 
%   the extrems are left intact.
%
%   Y = NANMOVING_AVERAGE(X,F,'interpNaN') Removes also the NaN element 
%   that has at least one NaN element at his surrounding.
%
%   Example:
%      x = 2*pi*linspace(-1,1); 
%      yn = cos(x) + 0.25 - 0.5*rand(size(x)); yn([5 30 40:43]) = NaN;
%      ys = nanmoving_average(yn,4,'interpnan');
%      plot(x,yn,x,ys), legend('noisy','smooth',4), axis tight
%
%   See also MOVING_AVERAGE, MOVING_AVERAGE2, NANMOVING_AVERAGE2

%   Written by
%   M.S. Carlos Adrián Vargas Aguilera
%   Physical Oceanography PhD candidate
%   CICESE 
%   Mexico, october 2006
%
%   nubeobscura@hotmail.com

suavenans = 0;
if (nargin == 3) && strcmpi(varargin{1}(1),'i')
 suavenans = 1;
end

% Moving average, except the ends:
Y = nanboxcar_win(X,F);

% Smooths the ends:
Wwidth = 2*F+1;
N = length(X);
Yini = nancummean(X(1:Wwidth-2));  
Y(1:F) = Yini(1:2:end);
Yfin = nancummean(X(N:-1:N-Wwidth+3));
Y(N-F+1:N) = Yfin(end:-2:1);

% Smooths NaN's?
if ~suavenans
 Y(isnan(X)) = NaN;
end


function Y = nancummean(X)
% Performs the acumulative mean of the vector ignoring NaN's.
%
% nubeobscura@hotmail.com

Y = X;
for n = 1:length(X)
 Y(n) = nan_mean( X(n:-1:1) );   
end

function Y = nan_mean(X)
% Performs the mean of a vector ignoring nans.
%
% nubeobscura@hotmail.com

X(isnan(X)) = [];
if isempty(X)
 Y = NaN;
else
 Y = mean(X);
end

function Y = nanboxcar_win(X,F)
% Boxcar window of width 2F+1 via recursive moving average and ignoring
% NaN's.
%
% nubeobscura@hotmail.com

Y = X;
if F == 0, return, end
N = length(X);
for n = F+1:N-F
 Y(n) = nan_mean( X(n-F:n+F) );
end


% Carlos Adrián Vargas Aguilera. nubeobscura@hotmail.com