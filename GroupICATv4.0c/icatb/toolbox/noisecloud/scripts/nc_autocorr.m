function varargout = nc_autocorr(Series , nLags , Q , nSTDs)
%AUTOCORR Compute or plot sample auto-correlation function.
%   Compute or plot the sample auto-correlation function (ACF) of a univariate,
%   stochastic time series. When called with no output arguments, AUTOCORR
%   displays the ACF sequence with confidence bounds.
%
%   [ACF, Lags, Bounds] = autocorr(Series)
%   [ACF, Lags, Bounds] = autocorr(Series , nLags , M , nSTDs)
%
%   Optional Inputs: nLags , M , nSTDs
%
% Inputs:
%   Series - Vector of observations of a univariate time series for which the
%     sample ACF is computed or plotted. The last row of Series contains the
%     most recent observation of the stochastic sequence.
%
% Optional Inputs:
%   nLags - Positive, scalar integer indicating the number of lags of the ACF
%     to compute. If empty or missing, the default is to compute the ACF at
%     lags 0,1,2, ... T = minimum[20 , length(Series)-1]. Since an ACF is
%     symmetric about zero lag, negative lags are ignored.
%
%   M - Non-negative integer scalar indicating the number of lags beyond which
%     the theoretical ACF is deemed to have died out. Under the hypothesis that
%     the underlying Series is really an MA(M) process, the large-lag standard
%     error is computed (via Bartlett's approximation) for lags > M as an
%     indication of whether the ACF is effectively zero beyond lag M. On the
%     assumption that the ACF is zero beyond lag M, Bartlett's approximation
%     is used to compute the standard deviation of the ACF for lags > M. If M
%     is empty or missing, the default is M = 0, in which case Series is
%     assumed to be Gaussian white noise. If Series is a Gaussian white noise
%     process of length N, the standard error will be approximately 1/sqrt(N).
%     M must be less than nLags.
%
%   nSTDs - Positive scalar indicating the number of standard deviations of the
%     sample ACF estimation error to compute assuming the theoretical ACF of
%     Series is zero beyond lag M. When M = 0 and Series is a Gaussian white
%     noise process of length N, specifying nSTDs will result in confidence
%     bounds at +/-(nSTDs/sqrt(N)). If empty or missing, default is nSTDs = 2
%     (i.e., approximate 95% confidence interval).
%
% Outputs:
%   ACF - Sample auto-correlation function of Series. ACF is a vector of
%     length nLags + 1 corresponding to lags 0,1,2,...,nLags. The first
%     element of ACF is unity (i.e., ACF(1) = 1 = lag 0 correlation).
%
%   Lags - Vector of lags corresponding to ACF (0,1,2,...,nLags).
%
%   Bounds - Two element vector indicating the approximate upper and lower
%     confidence bounds assuming that Series is an MA(M) process. Note that
%     Bounds is approximate for lags > M only.

[rows , columns]  =  size(Series);

if (rows ~= 1) & (columns ~= 1)
   error('GARCH:autocorr:NonVectorInput' , ' Input ''Series'' must be a vector.');
end

rowSeries   =  size(Series,1) == 1;

Series      =  Series(:);       % Ensure a column vector
n           =  length(Series);  % Sample size.
defaultLags =  20;              % BJR recommend about 20 lags for ACFs.

%
% Ensure the number of lags, nLags, is a positive
% integer scalar and set default if necessary.
%

if (nargin >= 2) & ~isempty(nLags)
  if prod(size(nLags)) > 1
     error('GARCH:autocorr:NonScalarLags' , ' Number of lags ''nLags'' must be a scalar.');
  end
  if (round(nLags) ~= nLags) | (nLags <= 0)
     error('GARCH:autocorr:NonPositiveInteger' , ' Number of lags ''nLags'' must be a positive integer.');
  end
  if nLags > (n - 1)
     error('GARCH:autocorr:LagsTooLarge' , ' Number of lags ''nLags'' must not exceed ''Series'' length - 1.');
  end
else
  nLags  =  min(defaultLags , n - 1);
end

%
% Ensure the hypothesized number of lags, Q, is a non-negative integer
% scalar, and set default if necessary.
%
if (nargin >= 3) & ~isempty(Q)
  if prod(size(Q)) > 1
     error('GARCH:autocorr:NonScalarQ' , ' Number of lags ''Q'' must be a scalar.');
  end
  if (round(Q) ~= Q) | (Q < 0)
     error('GARCH:autocorr:NegativeInteger' , ' Number of lags ''Q'' must be a non-negative integer.');
  end
  if Q >= nLags
     error('GARCH:autocorr:QTooLarge' , ' ''Q'' must be less than ''nLags''.');
  end
else
  Q  =  0;     % Default is 0 (Gaussian white noise hypothisis).
end

%
% Ensure the number of standard deviations, nSTDs, is a positive
% scalar and set default if necessary.
%

if (nargin >= 4) & ~isempty(nSTDs)
  if prod(size(nSTDs)) > 1
     error('GARCH:autocorr:NonScalarSTDs' , ' Number of standard deviations ''nSTDs'' must be a scalar.');
  end
  if nSTDs < 0
     error('GARCH:autocorr:NegativeSTDs' , ' Number of standard deviations ''nSTDs'' must be non-negative.');
  end
else
  nSTDs =  2;     % Default is 2 standard errors (95% condfidence interval).
end

%
% Convolution, polynomial multiplication, and FIR digital filtering are
% all the same operation. The FILTER command could be used to compute
% the ACF (by computing the correlation by convolving the de-meaned
% Series with a flipped version of itself), but FFT-based computation
% is significantly faster for large data sets.
%
% The ACF computation is based on Box, Jenkins, Reinsel, pages 30-34, 188.
%

nFFT =  2^(nextpow2(length(Series)) + 1);
F    =  fft(Series-mean(Series) , nFFT);
F    =  F .* conj(F);
ACF  =  ifft(F);
ACF  =  ACF(1:(nLags + 1));         % Retain non-negative lags.
ACF  =  ACF ./ ACF(1);     % Normalize.
ACF  =  real(ACF);

%
% Compute approximate confidence bounds using the Box-Jenkins-Reinsel
% approach, equations 2.1.13 and 6.2.2, on pages 33 and 188, respectively.
%

sigmaQ  =  sqrt((1 + 2*(ACF(2:Q+1)'*ACF(2:Q+1)))/n);
bounds  =  sigmaQ * [nSTDs ; -nSTDs];
Lags    =  [0:nLags]';

if nargout == 0                     % Make plot if requested.

%
%  Plot the sample ACF.
%
  lineHandles  =  stem(Lags , ACF , 'filled' , 'r-o');
  set   (lineHandles(1) , 'MarkerSize' , 4)
  grid  ('on')
  xlabel('Lag')
  ylabel('Sample Autocorrelation')
  title ('Sample Autocorrelation Function (ACF)')
  hold  ('on')
%
%  Plot the confidence bounds under the hypothesis that the underlying
%  Series is really an MA(Q) process. Bartlett's approximation gives
%  an indication of whether the ACF is effectively zero beyond lag Q.
%  For this reason, the confidence bounds (horizontal lines) appear
%  over the ACF ONLY for lags GREATER than Q (i.e., Q+1, Q+2, ... nLags).
%  In other words, the confidence bounds enclose ONLY those lags for
%  which the null hypothesis is assumed to hold.
%

  plot([Q+0.5 Q+0.5 ; nLags nLags] , [bounds([1 1]) bounds([2 2])] , '-b');

  plot([0 nLags] , [0 0] , '-k');
  hold('off')
  a  =  axis;
  axis([a(1:3) 1]);

else

%
%  Re-format outputs for compatibility with the SERIES input. When SERIES is
%  input as a row vector, then pass the outputs as a row vectors; when SERIES
%  is a column vector, then pass the outputs as a column vectors.
%
  if rowSeries
     ACF     =  ACF.';
     Lags    =  Lags.';
     bounds  =  bounds.';
  end

  varargout  =  {ACF , Lags , bounds};

end