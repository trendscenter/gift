function [Out1, Out2, Out3] = icatb_fastICA(mixedsig, varargin)
%FASTICA - Fast Independent Component Analysis
%
% FastICA for Matlab 5.x
% Version 2.1, January 15 2001
% Copyright (c) Hugo Gävert, Jarmo Hurri, Jaakko Särel? and Aapo Hyvärinen.
%
% FASTICA(mixedsig) estimates the independent components from given
% multidimensional signals. Each row of matrix mixedsig is one
% observed signal.  FASTICA uses Hyvarinen's fixed-point algorithm,
% see http://www.cis.hut.fi/projects/ica/fastica/. Output from the
% function depends on the number output arguments:
%
% [icasig] = FASTICA (mixedsig); the rows of icasig contain the
% estimated independent components.
%
% [icasig, A, W] = FASTICA (mixedsig); outputs the estimated separating
% matrix W and the corresponding mixing matrix A.
%
% [A, W] = FASTICA (mixedsig); gives only the estimated mixing matrix
% A and the separating matrix W.
%
% Some optional arguments induce other output formats, see below.
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
% FASTICA can be called with numerous optional arguments. Optional
% arguments are given in parameter pairs, so that first argument is
% the name of the parameter and the next argument is the value for
% that parameter. Optional parameter pairs can be given in any order.
%
% OPTIONAL PARAMETERS:
%
% Parameter name        Values and description
%
%======================================================================
% --Basic parameters in fixed-point algorithm:
%
% 'approach'            (string) The decorrelation approach used. Can be
%                       symmetric ('symm'), i.e. estimate all the
%                       independent component in parallel, or
%                       deflation ('defl'), i.e. estimate independent
%                       component one-by-one like in projection pursuit.
%                       Default is 'defl'.
%
% 'numOfIC'             (integer) Number of independent components to
%                       be estimated. Default equals the dimension of data.
%
%======================================================================
% --Choosing the nonlinearity:
%
% 'g'                   (string) Chooses the nonlinearity g used in 
%                       the fixed-point algorithm. Possible values:
%
%                       Value of 'g':      Nonlinearity used:
%                       'pow3' (default)   g(u)=u^3
%                       'tanh'             g(u)=tanh(a1*u)
%                       'gauss             g(u)=u*exp(-a2*u^2/2)
%                       'skew'             g(u)=u^2
% 
% 'finetune'		(string) Chooses the nonlinearity g used when 
%                       fine-tuning. In addition to same values
%                       as for 'g', the possible value 'finetune' is:
%                       'off'              fine-tuning is disabled.
%
% 'a1'                  (number) Parameter a1 used when g='tanh'.
%                       Default is 1.
% 'a2'                  (number) Parameter a2 used when g='gaus'.
%                       Default is 1.
%
% 'mu'			(number) Step size. Default is 1.
%                       If the value of mu is other than 1, then the
%                       program will use the stabilized version of the
%                       algorithm (see also parameter 'stabilization').
%
%
% 'stabilization'       (string) Values 'on' or 'off'. Default 'off'. 
%                       This parameter controls wether the program uses
%                       the stabilized version of the algorithm or
%                       not. If the stabilization is on, then the value
%                       of mu can momentarily be halved if the program
%                       senses that the algorithm is stuck between two
%                       points (this is called a stroke). Also if there
%                       is no convergence before half of the maximum
%                       number of iterations has been reached then mu
%                       will be halved for the rest of the rounds.
% 
%======================================================================
% --Controlling convergence:
%
% 'epsilon'             (number) Stopping criterion. Default is 0.0001.
%
% 'maxNumIterations'    (integer) Maximum number of iterations.
%                       Default is 1000.
%
% 'maxFinetune'         (integer) Maximum number of iterations in 
%                       fine-tuning. Default 100.
%
% 'sampleSize'          (number) [0 - 1] Percentage of samples used in
%                       one iteration. Samples are chosen in random.
%                       Default is 1 (all samples).
%
% 'initGuess'           (matrix) Initial guess for A. Default is random.
%                       You can now do a "one more" like this: 
%                       [ica, A, W] = fastica(mix, 'numOfIC',3);
%                       [ica2, A2, W2] = fastica(mix, 'initGuess', A, 'numOfIC', 4);
%
%======================================================================
% --Graphics and text output:
%
% 'verbose'             (string) Either 'on' or 'off'. Default is
%                       'on': report progress of algorithm in text format.
%
% 'displayMode'         (string) Plot running estimates of independent
%                       components: 'signals', 'basis', 'filters' or
%                       'off'. Default is 'off'~.    ~modifed by Yiou Li
%                       1/22/04
%
% 'displayInterval'     Number of iterations between plots.
%                       Default is 1 (plot after every iteration).
%
%======================================================================
% --Controlling reduction of dimension and whitening:
%
% Reduction of dimension is controlled by 'firstEig' and 'lastEig', or
% alternatively by 'interactivePCA'. 
%
% 'firstEig'            (integer) This and 'lastEig' specify the range for
%                       eigenvalues that are retained, 'firstEig' is
%                       the index of largest eigenvalue to be
%                       retained. Default is 1.
%
% 'lastEig'             (integer) This is the index of the last (smallest)
%                       eigenvalue to be retained. Default equals the
%                       dimension of data.
%
% 'interactivePCA'      (string) Either 'on' or 'off'. When set 'on', the
%                       eigenvalues are shown to the user and the
%                       range can be specified interactively. Default
%                       is 'off'. Can also be set to 'gui'. Then the user
%                       can use the same GUI that's in FASTICAG.
%
% If you already know the eigenvalue decomposition of the covariance
% matrix, you can avoid computing it again by giving it with the
% following options:
%
% 'pcaE'                (matrix) Eigenvectors
% 'pcaD'                (matrix) Eigenvalues
%
% If you already know the whitened data, you can give it directly to
% the algorithm using the following options:
%
% 'whiteSig'            (matrix) Whitened signal
% 'whiteMat'            (matrix) Whitening matrix
% 'dewhiteMat'          (matrix) dewhitening matrix
%
% If values for all the 'whiteSig', 'whiteSig' and 'dewhiteMat' are
% suplied, they will be used in computing the ICA. PCA and whitening
% are not performed. Though 'mixedsig' is not used in the main
% algorithm it still must be entered - some values are still
% calculated from it.
%
% Performing preprocessing only is possible by the option:
%
% 'only'                (string) Compute only PCA i.e. reduction of
%                       dimension ('pca') or only PCA plus whitening
%                       ('white'). Default is 'all': do ICA estimation
%                       as well.  This option changes the output
%                       format accordingly. For example: 
%
%                       [whitesig, WM, DWM] = FASTICA(mixedsig, 
%                       'only', 'white') 
%                       returns the whitened signals, the whitening matrix
%                       (WM) and the dewhitening matrix (DWM). (See also
%                       WHITENV.) In FastICA the whitening matrix performs
%                       whitening and the reduction of dimension. Dewhitening
%                       matrix is the pseudoinverse of whitening matrix.
%                        
%                       [E, D] = FASTICA(mixedsig, 'only', 'pca') 
%                       returns the eigenvector (E) and diagonal 
%                       eigenvalue (D) matrices  containing the 
%                       selected subspaces. 
%
%======================================================================
% EXAMPLES
%
%       [icasig] = FASTICA (mixedsig, 'approach', 'symm', 'g', 'tanh');
%               Do ICA with tanh nonlinearity and in parallel (like
%               maximum likelihood estimation for supergaussian data).
%
%       [icasig] = FASTICA (mixedsig, 'lastEig', 10, 'numOfIC', 3);
%               Reduce dimension to 10, and estimate only 3
%               independent components.
%
%       [icasig] = FASTICA (mixedsig, 'verbose', 'off', 'displayMode', 'off');
%               Don't output convergence reports and don't plot
%               independent components.
%
%
% A graphical user interface for FASTICA can be launched by the
% command FASTICAG
%
%   See also FASTICAG

% 15.1.2001
% Hugo Gävert


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the mean and check the data

[mixedsig, mixedmean] = remmean(mixedsig);

[Dim, NumOfSampl] = size(mixedsig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values for optional parameters

% All
verbose           = 'on';

% Default values for 'pcamat' parameters
firstEig          = 1;
lastEig           = Dim;
interactivePCA    = 'off';

% Default values for 'fpica' parameters
approach          = 'defl';
numOfIC           = Dim;
g                 = 'pow3';
finetune          = 'off';
a1                = 1;
a2                = 1;
myy               = 1;
stabilization     = 'off';
epsilon           = 0.0001;
maxNumIterations  = 1000;
maxFinetune       = 5;
initState         = 'rand';
guess             = 0;
sampleSize        = 1;
displayMode       = 'off';
displayInterval   = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters for fastICA - i.e. this file

b_verbose = 1;
jumpPCA = 0;
jumpWhitening = 0;
only = 3;
userNumOfIC = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the optional parameters

if (rem(length(varargin),2)==1)
  error('Optional parameters should always go by pairs');
else
  for i=1:2:(length(varargin)-1)
    % change the value of parameter
    switch varargin{i}
     case 'stabilization'
      stabilization = varargin{i+1};
     case 'maxFinetune'
      maxFinetune = varargin{i+1};
     case 'sampleSize'
      sampleSize = varargin{i+1};
     case 'verbose'
      verbose = varargin{i+1};
      % silence this program also
      if strcmp (verbose, 'off'), b_verbose = 0; end
     case 'firstEig'
      firstEig = varargin{i+1};
     case 'lastEig'
      lastEig = varargin{i+1};
     case 'interactivePCA'
      interactivePCA = varargin{i+1};
     case 'approach'
      approach = varargin{i+1};
     case 'numOfIC'
      numOfIC = varargin{i+1};
      % User has suplied new value for numOfIC.
      % We'll use this information later on...
      userNumOfIC = 1;
     case 'g'
      g = varargin{i+1};
     case 'finetune'
      finetune = varargin{i+1};
     case 'a1'
      a1 = varargin{i+1};
     case 'a2'
      a2 = varargin{i+1};
     case {'mu', 'myy'}
      myy = varargin{i+1};
     case 'epsilon'
      epsilon = varargin{i+1};
     case 'maxNumIterations'
      maxNumIterations = varargin{i+1};
     case 'initGuess'
      % no use setting 'guess' if the 'initState' is not set
      initState = 'guess';
      guess = varargin{i+1};
     case 'displayMode'
      displayMode = varargin{i+1};
     case 'displayInterval'
      displayInterval = varargin{i+1};
     case 'pcaE'
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      E = varargin{i+1};
     case 'pcaD'
      % calculate if there are enought parameters to skip PCA
      jumpPCA = jumpPCA + 1;
      D = varargin{i+1};
     case 'whiteSig'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whitesig = varargin{i+1};
     case 'whiteMat'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      whiteningMatrix = varargin{i+1};
     case 'dewhiteMat'
      % calculate if there are enought parameters to skip PCA and whitening
      jumpWhitening = jumpWhitening + 1;
      dewhiteningMatrix = varargin{i+1};
     case 'only'
      % if the user only wants to calculate PCA or...
      switch varargin{i+1}
       case 'pca'
	only = 1;
       case 'white'
	only = 2;
       case 'all'
	only = 3;
      end
      
     otherwise
      % Hmmm, something wrong with the parameter string
      error(['Unrecognized parameter: ''' varargin{i} '''']);
    end;
  end;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print information about data
if b_verbose
  fprintf('Number of signals: %d\n', Dim);
  fprintf('Number of samples: %d\n', NumOfSampl);
end

% Check if the data has been entered the wrong way,
% but warn only... it may be on purpose

if Dim > NumOfSampl
  if b_verbose
    fprintf('Warning: ');
    fprintf('The signal matrix may be oriented in the wrong way.\n');
    fprintf('In that case transpose the matrix.\n\n');
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating PCA

% We need the results of PCA for whitening, but if we don't
% need to do whitening... then we dont need PCA...
if jumpWhitening == 3
  if b_verbose,
    fprintf ('Whitened signal and corresponding matrises suplied.\n');
    fprintf ('PCA calculations not needed.\n');
  end;
else
  
  % OK, so first we need to calculate PCA
  % Check to see if we already have the PCA data
  if jumpPCA == 2,
    if b_verbose,
      fprintf ('Values for PCA calculations suplied.\n');
      fprintf ('PCA calculations not needed.\n');
    end;
  else
    % display notice if the user entered one, but not both, of E and D.
    if (jumpPCA > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump PCA:\n');
      fprintf ('''pcaE'', ''pcaD''.\n');
    end;
    
    % Calculate PCA
    [E, D]=pcamat(mixedsig, firstEig, lastEig, interactivePCA, verbose);
  end
end

% skip the rest if user only wanted PCA
if only > 1
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Whitening the data
  
  % Check to see if the whitening is needed...
  if jumpWhitening == 3,
    if b_verbose,
      fprintf ('Whitening not needed.\n');
    end;
  else
    
    % Whitening is needed
    % display notice if the user entered some of the whitening info, but not all.
    if (jumpWhitening > 0) & (b_verbose),
      fprintf ('You must suply all of these in order to jump whitening:\n');
      fprintf ('''whiteSig'', ''whiteMat'', ''dewhiteMat''.\n');
    end;
    
    % Calculate the whitening
    [whitesig, whiteningMatrix, dewhiteningMatrix] = whitenv ...
						     (mixedsig, E, D, verbose);
  end
  
end % if only > 1

% skip the rest if user only wanted PCA and whitening
if only > 2
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculating the ICA
  
  % Check some parameters
  % The dimension of the data may have been reduced during PCA calculations.
  % The original dimension is calculated from the data by default, and the
  % number of IC is by default set to equal that dimension.
  
  Dim = size(whitesig, 1);
  
  % The number of IC's must be less or equal to the dimension of data
  if numOfIC > Dim
    numOfIC = Dim;
    % Show warning only if verbose = 'on' and user suplied a value for 'numOfIC'
    if (b_verbose & userNumOfIC)
      fprintf('Warning: estimating only %d independent components\n', numOfIC);
      fprintf('(Can''t estimate more independent components than dimension of data)\n');
    end
  end
  
  % Calculate the ICA with fixed point algorithm.
  [A, W] = fpica (whitesig,  whiteningMatrix, dewhiteningMatrix, approach, ...
		  numOfIC, g, finetune, a1, a2, myy, stabilization, epsilon, ...
		  maxNumIterations, maxFinetune, initState, guess, sampleSize, ...
		  displayMode, displayInterval, verbose);
  
  % Check for valid return
  if ~isempty(W)
    % Add the mean back in.
    if b_verbose
      fprintf('Adding the mean back to the data.\n');
    end
    icasig = W * mixedsig + (W * mixedmean) * ones(1, NumOfSampl);
    %icasig = W * mixedsig;
    if b_verbose & ...
	  (max(abs(W * mixedmean)) > 1e-9) & ...
	  (strcmp(displayMode,'signals') | strcmp(displayMode,'on'))
      fprintf('Note that the plots don''t have the mean added.\n');
    end
  else
    icasig = [];
  end

end % if only > 2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The output depends on the number of output parameters
% and the 'only' parameter.

if only == 1    % only PCA
  Out1 = E;
  Out2 = D;
elseif only == 2  % only PCA & whitening
  if nargout == 2
    Out1 = whiteningMatrix;
    Out2 = dewhiteningMatrix;
  else
    Out1 = whitesig;
    Out2 = whiteningMatrix;
    Out3 = dewhiteningMatrix;
  end
else      % ICA
  if nargout == 2
    Out1 = A;
    Out2 = W;
  else
    Out1 = icasig;
    Out2 = A;
    Out3 = W;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A, W] = fpica(X, whiteningMatrix, dewhiteningMatrix, approach, ...
			numOfIC, g, finetune, a1, a2, myy, stabilization, ...
			epsilon, maxNumIterations, maxFinetune, initState, ...
			guess, sampleSize, displayMode, displayInterval, ...
			s_verbose);
%FPICA - Fixed point ICA. Main algorithm of FASTICA.
%
% [A, W] = fpica(whitesig, whiteningMatrix, dewhiteningMatrix, approach,
%        numOfIC, g, finetune, a1, a2, mu, stabilization, epsilon, 
%        maxNumIterations, maxFinetune, initState, guess, sampleSize,
%        displayMode, displayInterval, verbose);
% 
% Perform independent component analysis using Hyvarinen's fixed point
% algorithm. Outputs an estimate of the mixing matrix A and its inverse W.
%
% whitesig                              :the whitened data as row vectors
% whiteningMatrix                       :can be obtained with function whitenv
% dewhiteningMatrix                     :can be obtained with function whitenv
% approach      [ 'symm' | 'defl' ]     :the approach used (deflation or symmetric)
% numOfIC       [ 0 - Dim of whitesig ] :number of independent components estimated
% g             [ 'pow3' | 'tanh' |     :the nonlinearity used
%                 'gaus' | 'skew' ]     
% finetune      [same as g + 'off']     :the nonlinearity used in finetuning.
% a1                                    :parameter for tuning 'tanh'
% a2                                    :parameter for tuning 'gaus'
% mu                                    :step size in stabilized algorithm
% stabilization [ 'on' | 'off' ]        :if mu < 1 then automatically on
% epsilon                               :stopping criterion
% maxNumIterations                      :maximum number of iterations 
% maxFinetune                           :maximum number of iteretions for finetuning
% initState     [ 'rand' | 'guess' ]    :initial guess or random initial state. See below
% guess                                 :initial guess for A. Ignored if initState = 'rand'
% sampleSize    [ 0 - 1 ]               :percentage of the samples used in one iteration
% displayMode   [ 'signals' | 'basis' | :plot running estimate
%                 'filters' | 'off' ]
% displayInterval                       :number of iterations we take between plots
% verbose       [ 'on' | 'off' ]        :report progress in text format
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%       [nv, wm, dwm] = whitenv(vectors, E, D);
%       [A, W] = fpica(nv, wm, dwm);
%
%
% This function is needed by FASTICA and FASTICAG
%
%   See also FASTICA, FASTICAG, WHITENV, PCAMAT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Global variable for stopping the ICA calculations from the GUI
global g_FastICA_interrupt;
if isempty(g_FastICA_interrupt)
  clear global g_FastICA_interrupt;
  interruptible = 0;
else
  interruptible = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values

if nargin < 3, error('Not enough arguments!'); end
[vectorSize, numSamples] = size(X);
if nargin < 20, s_verbose = 'on'; end
if nargin < 19, displayInterval = 1; end
if nargin < 18, displayMode = 'on'; end
if nargin < 17, sampleSize = 1; end
if nargin < 16, guess = 1; end
if nargin < 15, initState = 'rand'; end
if nargin < 14, maxFinetune = 100; end
if nargin < 13, maxNumIterations = 1000; end
if nargin < 12, epsilon = 0.0001; end
if nargin < 11, stabilization = 'on'; end
if nargin < 10, myy = 1; end
if nargin < 9, a2 = 1; end
if nargin < 8, a1 = 1; end
if nargin < 7, finetune = 'off'; end
if nargin < 6, g = 'pow3'; end
if nargin < 5, numOfIC = vectorSize; end     % vectorSize = Dim
if nargin < 4, approach = 'defl'; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the data

if ~isreal(X)
  error('Input has an imaginary part.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for verbose

switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for approach

switch lower(approach)
 case 'symm'
  approachMode = 1;
 case 'defl'
  approachMode = 2;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''approach''\n', approach));
end
if b_verbose, fprintf('Used approach [ %s ].\n', approach); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for numOfIC

if vectorSize < numOfIC
  error('Must have numOfIC <= Dimension!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the sampleSize
if sampleSize > 1
  sampleSize = 1;
  if b_verbose
    fprintf('Warning: Setting ''sampleSize'' to 1.\n');
  end  
elseif sampleSize < 1
  if (sampleSize * numSamples) < 1000
    sampleSize = min(1000/numSamples, 1);
    if b_verbose
      fprintf('Warning: Setting ''sampleSize'' to %0.3f (%d samples).\n', ...
	      sampleSize, floor(sampleSize * numSamples));
    end  
  end
end
if b_verbose
  if  b_verbose & (sampleSize < 1)
    fprintf('Using about %0.0f%% of the samples in random order in every step.\n',sampleSize*100);
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for nonlinearity.

switch lower(g)
 case 'pow3'
  gOrig = 10;
 case 'tanh'
  gOrig = 20;
 case {'gaus', 'gauss'}
  gOrig = 30;
 case 'skew'
  gOrig = 40;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''g''\n', g));
end
if sampleSize ~= 1
  gOrig = gOrig + 2;
end
if myy ~= 1
  gOrig = gOrig + 1;
end

if b_verbose,
  fprintf('Used nonlinearity [ %s ].\n', g);
end

finetuningEnabled = 1;
switch lower(finetune)
 case 'pow3'
  gFine = 10 + 1;
 case 'tanh'
  gFine = 20 + 1;
 case {'gaus', 'gauss'}
  gFine = 30 + 1;
 case 'skew'
  gFine = 40 + 1;
 case 'off'
  if myy ~= 1
    gFine = gOrig;
  else 
    gFine = gOrig + 1;
  end
  finetuningEnabled = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''finetune''\n', ...
		finetune));
end

if b_verbose & finetuningEnabled
  fprintf('Finetuning enabled (nonlinearity: [ %s ]).\n', finetune);
end

switch lower(stabilization)
 case 'on'
  stabilizationEnabled = 1;
 case 'off'
  if myy ~= 1
    stabilizationEnabled = 1;
  else
    stabilizationEnabled = 0;
  end
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''stabilization''\n', ...
		stabilization)); 
end

if b_verbose & stabilizationEnabled
  fprintf('Using stabilized algorithm.\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some other parameters
myyOrig = myy;
% When we start fine-tuning we'll set myy = myyK * myy
myyK = 0.01;
% How many times do we try for convergence until we give up.
failureLimit = 5;


usedNlinearity = gOrig;
stroke = 0;
notFine = 1;
long = 0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for initial state.

switch lower(initState)
 case 'rand'
  initialStateMode = 0;
 case 'guess'
  if size(guess,1) ~= size(whiteningMatrix,2)
    initialStateMode = 0;
    if b_verbose
      fprintf('Warning: size of initial guess is incorrect. Using random initial guess.\n');
    end
  else
    initialStateMode = 1;
    if size(guess,2) < numOfIC
      if b_verbose
	fprintf('Warning: initial guess only for first %d components. Using random initial guess for others.\n', size(guess,2)); 
      end
      guess(:, size(guess, 2) + 1:numOfIC) = ...
					     rand(vectorSize,numOfIC-size(guess,2))-.5;
    elseif size(guess,2)>numOfIC
      guess=guess(:,1:numOfIC);
      fprintf('Warning: Initial guess too large. The excess column are dropped.\n');
    end
    if b_verbose, fprintf('Using initial guess.\n'); end
  end
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''initState''\n', initState));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Checking the value for display mode.

switch lower(displayMode)
 case {'off', 'none'}
  usedDisplay = 0;
 case {'on', 'signals'}
  usedDisplay = 1;
  if (b_verbose & (numSamples > 10000))
    fprintf('Warning: Data vectors are very long. Plotting may take long time.\n');
  end
  if (b_verbose & (numOfIC > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
 case 'basis'
  usedDisplay = 2;
  if (b_verbose & (numOfIC > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
 case 'filters'
  usedDisplay = 3;
  if (b_verbose & (vectorSize > 25))
    fprintf('Warning: There are too many signals to plot. Plot may not look good.\n');
  end
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''displayMode''\n', displayMode));
end

% The displayInterval can't be less than 1...
if displayInterval < 1
  displayInterval = 1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if b_verbose, fprintf('Starting ICA calculation...\n'); end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SYMMETRIC APPROACH
if approachMode == 1,

  % set some parameters more...
  usedNlinearity = gOrig;
  stroke = 0;
  notFine = 1;
  long = 0;
  
  A = zeros(vectorSize, numOfIC);  % Dewhitened basis vectors.
  if initialStateMode == 0
    % Take random orthonormal initial vectors.
    B = orth(rand(vectorSize, numOfIC) - .5);
  elseif initialStateMode == 1
    % Use the given initial vector as the initial state
    B = whiteningMatrix * guess;
  end
  
  BOld = zeros(size(B));
  BOld2 = zeros(size(B));
  
  % This is the actual fixed-point iteration loop.
  for round = 1:maxNumIterations + 1,
    if round == maxNumIterations + 1,
      if b_verbose 
        fprintf('No convergence after %d steps\n', maxNumIterations);
        fprintf('Note that the plots are probably wrong.\n');
      end
      A=[];
      W=[];
      return;
    end
    
    if (interruptible & g_FastICA_interrupt)
      if b_verbose 
        fprintf('\n\nCalculation interrupted by the user\n');
      end
      if ~isempty(B)
	W = B' * whiteningMatrix;
	A = dewhiteningMatrix * B;
      else
	W = [];
	A = [];
      end
      return;
    end
    
    
    % Symmetric orthogonalization.
    B = B * real(inv(B' * B)^(1/2));
    
    % Test for termination condition. Note that we consider opposite
    % directions here as well.
    minAbsCos = min(abs(diag(B' * BOld)));
    minAbsCos2 = min(abs(diag(B' * BOld2)));
    
    if (1 - minAbsCos < epsilon)
      if finetuningEnabled & notFine
        if b_verbose, fprintf('Initial convergence, fine-tuning: \n'); end;
        notFine = 0;
        usedNlinearity = gFine;
        myy = myyK * myyOrig;
        BOld = zeros(size(B));
        BOld2 = zeros(size(B));
	
      else
        if b_verbose, fprintf('Convergence after %d steps\n', round); end
	
        % Calculate the de-whitened vectors.
        A = dewhiteningMatrix * B;
        break;
      end
    elseif stabilizationEnabled
      if (~stroke) & (1 - minAbsCos2 < epsilon)
	if b_verbose, fprintf('Stroke!\n'); end;
	stroke = myy;
	myy = .5*myy;
	if mod(usedNlinearity,2) == 0
	  usedNlinearity = usedNlinearity + 1;
	end
      elseif stroke
	myy = stroke;
	stroke = 0;
	if (myy == 1) & (mod(usedNlinearity,2) ~= 0)
	  usedNlinearity = usedNlinearity - 1;
	end
      elseif (~long) & (round>maxNumIterations/2)
	if b_verbose, fprintf('Taking long (reducing step size)\n'); end;
	long = 1;
	myy = .5*myy;
	if mod(usedNlinearity,2) == 0
	  usedNlinearity = usedNlinearity + 1;
	end
      end
    end
    
    BOld2 = BOld;
    BOld = B;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show the progress...
    if b_verbose
      if round == 1
        fprintf('Step no. %d\n', round);
      else
        fprintf('Step no. %d, change in value of estimate: %.3f \n', round, 1 - minAbsCos);
      end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Also plot the current state...
    switch usedDisplay
     case 1
      if rem(round, displayInterval) == 0,
	% There was and may still be other displaymodes...
	% 1D signals
	icaplot('dispsig',(X'*B)');
      end
     case 2
      if rem(round, displayInterval) == 0,
	% ... and now there are :-)
	% 1D basis
	A = dewhiteningMatrix * B;
	icaplot('dispsig',A');
      end
     case 3
      if rem(round, displayInterval) == 0,
	% ... and now there are :-)
	% 1D filters
	W = B' * whiteningMatrix;
	icaplot('dispsig',W);
      end
     otherwise
    end
    
    drawnow;
    
    switch usedNlinearity
      % pow3
     case 10
      B = (X * (( X' * B) .^ 3)) / numSamples - 3 * B;
     case 11
      % optimoitu - epsilonin kokoisia eroja
      % täm?on optimoitu koodi, katso vanha koodi esim.
      % aikaisemmista versioista kuten 2.0 beta3
      Y = X' * B;
      Gpow3 = Y .^ 3;
      Beta = sum(Y .* Gpow3);
      D = diag(1 ./ (Beta - 3 * numSamples));
      B = B + myy * B * (Y' * Gpow3 - diag(Beta)) * D;
     case 12
      Xsub=X(:, getSamples(numSamples, sampleSize));
      B = (Xsub * (( Xsub' * B) .^ 3)) / size(Xsub,2) - 3 * B;
     case 13
      % Optimoitu
      Ysub=X(:, getSamples(numSamples, sampleSize))' * B;
      Gpow3 = Ysub .^ 3;
      Beta = sum(Ysub .* Gpow3);
      D = diag(1 ./ (Beta - 3 * size(Ysub', 2)));
      B = B + myy * B * (Ysub' * Gpow3 - diag(Beta)) * D;
      
      % tanh
     case 20
      hypTan = tanh(a1 * X' * B);
      B = X * hypTan / numSamples - ...
	  ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / numSamples * ...
	  a1;
     case 21
      % optimoitu - epsilonin kokoisia 
      Y = X' * B;
      hypTan = tanh(a1 * Y);
      Beta = sum(Y .* hypTan);
      D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2)));
      B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;
     case 22
      Xsub=X(:, getSamples(numSamples, sampleSize));
      hypTan = tanh(a1 * Xsub' * B);
      B = Xsub * hypTan / size(Xsub, 2) - ...
	  ones(size(B,1),1) * sum(1 - hypTan .^ 2) .* B / size(Xsub, 2) * a1;
     case 23
      % Optimoitu
      Y = X(:, getSamples(numSamples, sampleSize))' * B;
      hypTan = tanh(a1 * Y);
      Beta = sum(Y .* hypTan);
      D = diag(1 ./ (Beta - a1 * sum(1 - hypTan .^ 2)));
      B = B + myy * B * (Y' * hypTan - diag(Beta)) * D;
      
      % gauss
     case 30
      U = X' * B;
      Usquared=U .^ 2;
      ex = exp(-a2 * Usquared / 2);
      gauss =  U .* ex;
      dGauss = (1 - a2 * Usquared) .*ex;
      B = X * gauss / numSamples - ...
	  ones(size(B,1),1) * sum(dGauss)...
	  .* B / numSamples ;
     case 31
      % optimoitu
      Y = X' * B;
      ex = exp(-a2 * (Y .^ 2) / 2);
      gauss = Y .* ex;
      Beta = sum(Y .* gauss);
      D = diag(1 ./ (Beta - sum((1 - a2 * (Y .^ 2)) .* ex)));
      B = B + myy * B * (Y' * gauss - diag(Beta)) * D;
     case 32
      Xsub=X(:, getSamples(numSamples, sampleSize));
      U = Xsub' * B;
      Usquared=U .^ 2;
      ex = exp(-a2 * Usquared / 2);
      gauss =  U .* ex;
      dGauss = (1 - a2 * Usquared) .*ex;
      B = Xsub * gauss / size(Xsub,2) - ...
	  ones(size(B,1),1) * sum(dGauss)...
	  .* B / size(Xsub,2) ;
     case 33
      % Optimoitu
      Y = X(:, getSamples(numSamples, sampleSize))' * B;
      ex = exp(-a2 * (Y .^ 2) / 2);
      gauss = Y .* ex;
      Beta = sum(Y .* gauss);
      D = diag(1 ./ (Beta - sum((1 - a2 * (Y .^ 2)) .* ex)));
      B = B + myy * B * (Y' * gauss - diag(Beta)) * D;
      
      % skew
     case 40
      B = (X * ((X' * B) .^ 2)) / numSamples;
     case 41
      % Optimoitu
      Y = X' * B;
      Gskew = Y .^ 2;
      Beta = sum(Y .* Gskew);
      D = diag(1 ./ (Beta));
      B = B + myy * B * (Y' * Gskew - diag(Beta)) * D;
     case 42
      Xsub=X(:, getSamples(numSamples, sampleSize));
      B = (Xsub * ((Xsub' * B) .^ 2)) / size(Xsub,2);
     case 43
      % Uusi optimoitu
      Y = X(:, getSamples(numSamples, sampleSize))' * B;
      Gskew = Y .^ 2;
      Beta = sum(Y .* Gskew);
      D = diag(1 ./ (Beta));
      B = B + myy * B * (Y' * Gskew - diag(Beta)) * D;

     otherwise
      error('Code for desired nonlinearity not found!');
    end
  end

  
  % Calculate ICA filters.
  W = B' * whiteningMatrix;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Also plot the last one...
  switch usedDisplay
   case 1 
    % There was and may still be other displaymodes...
    % 1D signals
    icaplot('dispsig',(X'*B)');
    drawnow;
   case 2
    % ... and now there are :-)
    % 1D basis
    icaplot('dispsig',A');
    drawnow;
   case 3
    % ... and now there are :-)
    % 1D filters
    icaplot('dispsig',W);
    drawnow;
   otherwise
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFLATION APPROACH
if approachMode == 2
  
  B = zeros(vectorSize);
  
  % The search for a basis vector is repeated numOfIC times.
  round = 1;
  
  numFailures = 0;
  
  while round <= numOfIC,
    myy = myyOrig;
    usedNlinearity = gOrig;
    stroke = 0;
    notFine = 1;
    long = 0;
    endFinetuning = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Show the progress...
    if b_verbose, fprintf('IC %d ', round); end
    
    % Take a random initial vector of lenght 1 and orthogonalize it
    % with respect to the other vectors.
    if initialStateMode == 0
      w = rand(vectorSize, 1) - .5;
    elseif initialStateMode == 1
      w=whiteningMatrix*guess(:,round);
    end
    w = w - B * B' * w;
    w = w / norm(w);
    
    wOld = zeros(size(w));
    wOld2 = zeros(size(w));
    
    % This is the actual fixed-point iteration loop.
    %    for i = 1 : maxNumIterations + 1
    i = 1;
    gabba = 1;
    while i <= maxNumIterations + gabba
      drawnow;
      if (interruptible & g_FastICA_interrupt)
        if b_verbose 
          fprintf('\n\nCalculation interrupted by the user\n');
        end
        return;
      end
      
      % Project the vector into the space orthogonal to the space
      % spanned by the earlier found basis vectors. Note that we can do
      % the projection with matrix B, since the zero entries do not
      % contribute to the projection.
      w = w - B * B' * w;
      w = w / norm(w);
      
      if notFine
	if i == maxNumIterations + 1
	  if b_verbose
	    fprintf('\nComponent number %d did not converge in %d iterations.\n', round, maxNumIterations);
	  end
	  round = round - 1;
	  numFailures = numFailures + 1;
	  if numFailures > failureLimit
	    if b_verbose
	      fprintf('Too many failures to converge (%d). Giving up.\n', numFailures);
	    end
	    if round == 0
	      A=[];
	      W=[];
	    end
	    return;
	  end
	  % numFailures > failurelimit
	  break;
	end
	% i == maxNumIterations + 1
      else
	% if notFine
	if i >= endFinetuning
	  wOld = w; % So the algorithm will stop on the next test...
	end
      end
      % if notFine
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Show the progress...
      if b_verbose, fprintf('.'); end;
      
      
      % Test for termination condition. Note that the algorithm has
      % converged if the direction of w and wOld is the same, this
      % is why we test the two cases.
      if norm(w - wOld) < epsilon | norm(w + wOld) < epsilon
        if finetuningEnabled & notFine
          if b_verbose, fprintf('Initial convergence, fine-tuning: '); end;
          notFine = 0;
	  gabba = maxFinetune;
          wOld = zeros(size(w));
          wOld2 = zeros(size(w));
          usedNlinearity = gFine;
          myy = myyK * myyOrig;
	  
	  endFinetuning = maxFinetune + i;
	  
        else
          numFailures = 0;
          % Save the vector
          B(:, round) = w;
	  
          % Calculate the de-whitened vector.
          A(:,round) = dewhiteningMatrix * w;
          % Calculate ICA filter.
          W(round,:) = w' * whiteningMatrix;
	  
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Show the progress...
          if b_verbose, fprintf('computed ( %d steps ) \n', i); end
	  
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          % Also plot the current state...
          switch usedDisplay
	   case 1
	    if rem(round, displayInterval) == 0,
	      % There was and may still be other displaymodes...   
	      % 1D signals
	      temp = X'*B;
	      icaplot('dispsig',temp(:,1:numOfIC)');
	      drawnow;
	    end
	   case 2
	    if rem(round, displayInterval) == 0,
	      % ... and now there are :-) 
	      % 1D basis
	      icaplot('dispsig',A');
	      drawnow;
	    end
	   case 3
	    if rem(round, displayInterval) == 0,
	      % ... and now there are :-) 
	      % 1D filters
	      icaplot('dispsig',W);
	      drawnow;
	    end
          end
	  % switch usedDisplay
	  break; % IC ready - next...
        end
	%if finetuningEnabled & notFine
      elseif stabilizationEnabled
	if (~stroke) & (norm(w - wOld2) < epsilon | norm(w + wOld2) < ...
			epsilon)
	  stroke = myy;
	  if b_verbose, fprintf('Stroke!'); end;
	  myy = .5*myy;
	  if mod(usedNlinearity,2) == 0
	    usedNlinearity = usedNlinearity + 1;
	  end
	elseif stroke
	  myy = stroke;
	  stroke = 0;
	  if (myy == 1) & (mod(usedNlinearity,2) ~= 0)
	    usedNlinearity = usedNlinearity - 1;
	  end
	elseif (notFine) & (~long) & (i > maxNumIterations / 2)
	  if b_verbose, fprintf('Taking long (reducing step size) '); end;
	  long = 1;
	  myy = .5*myy;
	  if mod(usedNlinearity,2) == 0
	    usedNlinearity = usedNlinearity + 1;
	  end
	end
      end
      
      wOld2 = wOld;
      wOld = w;
      
      switch usedNlinearity
	% pow3
       case 10
	w = (X * ((X' * w) .^ 3)) / numSamples - 3 * w;
       case 11
	EXGpow3 = (X * ((X' * w) .^ 3)) / numSamples;
	Beta = w' * EXGpow3;
	w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
       case 12
	Xsub=X(:,getSamples(numSamples, sampleSize));
	w = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2) - 3 * w;
       case 13
	Xsub=X(:,getSamples(numSamples, sampleSize));
	EXGpow3 = (Xsub * ((Xsub' * w) .^ 3)) / size(Xsub, 2);
	Beta = w' * EXGpow3;
	w = w - myy * (EXGpow3 - Beta * w) / (3 - Beta);
	% tanh
       case 20
	hypTan = tanh(a1 * X' * w);
	w = (X * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / numSamples;
       case 21
	hypTan = tanh(a1 * X' * w);
	Beta = w' * X * hypTan;
	w = w - myy * ((X * hypTan - Beta * w) / ...
		       (a1 * sum((1-hypTan .^2)') - Beta));
       case 22
	Xsub=X(:,getSamples(numSamples, sampleSize));
	hypTan = tanh(a1 * Xsub' * w);
	w = (Xsub * hypTan - a1 * sum(1 - hypTan .^ 2)' * w) / size(Xsub, 2);
       case 23
	Xsub=X(:,getSamples(numSamples, sampleSize));
	hypTan = tanh(a1 * Xsub' * w);
	Beta = w' * Xsub * hypTan;
	w = w - myy * ((Xsub * hypTan - Beta * w) / ...
		       (a1 * sum((1-hypTan .^2)') - Beta));
	% gauss
       case 30
	% This has been split for performance reasons.
	u = X' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	w = (X * gauss - sum(dGauss)' * w) / numSamples;
       case 31
	u = X' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	Beta = w' * X * gauss;
	w = w - myy * ((X * gauss - Beta * w) / ...
		       (sum(dGauss)' - Beta));
       case 32
	Xsub=X(:,getSamples(numSamples, sampleSize));
	u = Xsub' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	w = (Xsub * gauss - sum(dGauss)' * w) / size(Xsub, 2);
       case 33
	Xsub=X(:,getSamples(numSamples, sampleSize));
	u = Xsub' * w;
	u2=u.^2;
	ex=exp(-a2 * u2/2);
	gauss =  u.*ex;
	dGauss = (1 - a2 * u2) .*ex;
	Beta = w' * Xsub * gauss;
	w = w - myy * ((Xsub * gauss - Beta * w) / ...
		       (sum(dGauss)' - Beta));
	% skew
       case 40
	w = (X * ((X' * w) .^ 2)) / numSamples;
       case 41
	EXGskew = (X * ((X' * w) .^ 2)) / numSamples;
	Beta = w' * EXGskew;
	w = w - myy * (EXGskew - Beta*w)/(-Beta);
       case 42
	Xsub=X(:,getSamples(numSamples, sampleSize));
	w = (Xub * ((Xub' * w) .^ 2)) / size(Xsub, 2);
       case 43
	Xsub=X(:,getSamples(numSamples, sampleSize));
	EXGskew = (Xsub * ((Xsub' * w) .^ 2)) / size(Xsub, 2);
	Beta = w' * EXGskew;
	w = w - myy * (EXGskew - Beta*w)/(-Beta);
	
       otherwise
	error('Code for desired nonlinearity not found!');
      end
      
      % Normalize the new w.
      w = w / norm(w);
      i = i + 1;
    end
    round = round + 1;
  end
  if b_verbose, fprintf('Done.\n'); end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Also plot the ones that may not have been plotted.
  if (usedDisplay > 0) & (rem(round-1, displayInterval) ~= 0)
    switch usedDisplay
     case 1
      % There was and may still be other displaymodes...
      % 1D signals
      temp = X'*B;
      icaplot('dispsig',temp(:,1:numOfIC)');
      drawnow;
     case 2
      % ... and now there are :-)
      % 1D basis
      icaplot('dispsig',A');
      drawnow;
     case 3
      % ... and now there are :-)
      % 1D filters
      icaplot('dispsig',W);
      drawnow;
     otherwise
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the end let's check the data for some security
if ~isreal(A)
  if b_verbose, fprintf('Warning: removing the imaginary part from the result.\n'); end
  A = real(A);
  W = real(W);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subfunction
% Calculates tanh simplier and faster than Matlab tanh.
function y=tanh(x)
y = 1 - 2 ./ (exp(2 * x) + 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Samples = getSamples(max, percentage)
Samples = find(rand(1, max) < percentage);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [E, D] = pcamat(vectors, firstEig, lastEig, s_interactive, ...
    s_verbose);
%PCAMAT - Calculates the pca for data
%
% [E, D] = pcamat(vectors, firstEig, lastEig, ... 
%                 interactive, verbose);
%
% Calculates the PCA matrices for given data (row) vectors. Returns
% the eigenvector (E) and diagonal eigenvalue (D) matrices containing the
% selected subspaces. Dimensionality reduction is controlled with
% the parameters 'firstEig' and 'lastEig' - but it can also be done
% interactively by setting parameter 'interactive' to 'on' or 'gui'.
%
% ARGUMENTS
%
% vectors       Data in row vectors.
% firstEig      Index of the largest eigenvalue to keep.
%               Default is 1.
% lastEig       Index of the smallest eigenvalue to keep.
%               Default is equal to dimension of vectors.
% interactive   Specify eigenvalues to keep interactively. Note that if
%               you set 'interactive' to 'on' or 'gui' then the values
%               for 'firstEig' and 'lastEig' will be ignored, but they
%               still have to be entered. If the value is 'gui' then the
%               same graphical user interface as in FASTICAG will be
%               used. Default is 'off'.
% verbose       Default is 'on'.
%
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%
% Note 
%       The eigenvalues and eigenvectors returned by PCAMAT are not sorted.
%
% This function is needed by FASTICA and FASTICAG

% For historical reasons this version does not sort the eigenvalues or
% the eigen vectors in any ways. Therefore neither does the FASTICA or
% FASTICAG. Generally it seams that the components returned from
% whitening is almost in reversed order. (That means, they usually are,
% but sometime they are not - depends on the EIG-command of matlab.)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Default values:
if nargin < 5, s_verbose = 'on'; end
if nargin < 4, s_interactive = 'off'; end
if nargin < 3, lastEig = size(vectors, 1); end
if nargin < 2, firstEig = 1; end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check the optional parameters;
switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

switch lower(s_interactive)
 case 'on'
  b_interactive = 1;
 case 'off'
  b_interactive = 0;
 case 'gui'
  b_interactive = 2;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''interactive''\n', ...
		s_interactive));
end

oldDimension = size (vectors, 1);
if ~(b_interactive)
  if lastEig < 1 | lastEig > oldDimension
    error(sprintf('Illegal value [ %d ] for parameter: ''lastEig''\n', lastEig));
  end
  if firstEig < 1 | firstEig > lastEig
    error(sprintf('Illegal value [ %d ] for parameter: ''firstEig''\n', firstEig));
  end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate PCA

% Calculate the covariance matrix.
if b_verbose, fprintf ('Calculating covariance...\n'); end
covarianceMatrix = cov(vectors', 1);

maxLastEig = rank(covarianceMatrix, 1e-9);

% Calculate the eigenvalues and eigenvectors of covariance matrix.
[E, D] = eig(covarianceMatrix);

% Sort the eigenvalues - decending.
eigenvalues = flipud(sort(diag(D)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive part - command-line
if b_interactive == 1

  % Show the eigenvalues to the user
  hndl_win=figure;
  bar(eigenvalues);
  title('Eigenvalues');

  % ask the range from the user...
  % ... and keep on asking until the range is valid :-)
  areValuesOK=0;
  while areValuesOK == 0
    firstEig = input('The index of the largest eigenvalue to keep? (1) ');
    lastEig = input(['The index of the smallest eigenvalue to keep? (' ...
                    int2str(oldDimension) ') ']);
    % Check the new values...
    % if they are empty then use default values
    if isempty(firstEig), firstEig = 1;end
    if isempty(lastEig), lastEig = oldDimension;end
    % Check that the entered values are within the range
    areValuesOK = 1;
    if lastEig < 1 | lastEig > oldDimension
      fprintf('Illegal number for the last eigenvalue.\n');
      areValuesOK = 0;
    end
    if firstEig < 1 | firstEig > lastEig
      fprintf('Illegal number for the first eigenvalue.\n');
      areValuesOK = 0;
    end
  end
  % close the window
  close(hndl_win);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Interactive part - GUI
if b_interactive == 2

  % Show the eigenvalues to the user
  hndl_win = figure('Color',[0.8 0.8 0.8], ...
    'PaperType','a4letter', ...
    'Units', 'normalized', ...
    'Name', 'FastICA: Reduce dimension', ...
    'NumberTitle','off', ...
    'Tag', 'f_eig');
  h_frame = uicontrol('Parent', hndl_win, ...
    'BackgroundColor',[0.701961 0.701961 0.701961], ...
    'Units', 'normalized', ...
    'Position',[0.13 0.05 0.775 0.17], ...
    'Style','frame', ...
    'Tag','f_frame');

b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.142415 0.0949436 0.712077 0.108507], ...
	'String','Give the indices of the largest and smallest eigenvalues of the covariance matrix to be included in the reduced data.', ...
	'Style','text', ...
	'Tag','StaticText1');
e_first = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'f=round(str2num(get(gcbo, ''String'')));' ...
          'if (f < 1), f=1; end;' ...
          'l=str2num(get(findobj(''Tag'',''e_last''), ''String''));' ...
          'if (f > l), f=l; end;' ...
          'set(gcbo, ''String'', int2str(f));' ...
          ], ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'Position',[0.284831 0.0678168 0.12207 0.0542535], ...
	'Style','edit', ...
        'String', '1', ...
	'Tag','e_first');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.142415 0.0678168 0.12207 0.0542535], ...
	'String','Range from', ...
	'Style','text', ...
	'Tag','StaticText2');
e_last = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'l=round(str2num(get(gcbo, ''String'')));' ...
          'lmax = get(gcbo, ''UserData'');' ...
          'if (l > lmax), l=lmax; fprintf([''The selected value was too large, or the selected eigenvalues were close to zero\n'']); end;' ...
          'f=str2num(get(findobj(''Tag'',''e_first''), ''String''));' ...
          'if (l < f), l=f; end;' ...
          'set(gcbo, ''String'', int2str(l));' ...
          ], ...
	'BackgroundColor',[1 1 1], ...
	'HorizontalAlignment','right', ...
	'Position',[0.467936 0.0678168 0.12207 0.0542535], ...
	'Style','edit', ...
        'String', int2str(maxLastEig), ...
        'UserData', maxLastEig, ...
	'Tag','e_last');
% in the first version oldDimension was used instead of 
% maxLastEig, but since the program would automatically
% drop the eigenvalues afte maxLastEig...
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'BackgroundColor',[0.701961 0.701961 0.701961], ...
	'HorizontalAlignment','left', ...
	'Position',[0.427246 0.0678168 0.0406901 0.0542535], ...
	'String','to', ...
	'Style','text', ...
	'Tag','StaticText3');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback','uiresume(gcbf)', ...
	'Position',[0.630697 0.0678168 0.12207 0.0542535], ...
	'String','OK', ...
	'Tag','Pushbutton1');
b = uicontrol('Parent',hndl_win, ...
	'Units','normalized', ...
	'Callback',[ ...
          'gui_help(''pcamat'');' ...
          ], ...
	'Position',[0.767008 0.0678168 0.12207 0.0542535], ...
	'String','Help', ...
	'Tag','Pushbutton2');

  h_axes = axes('Position' ,[0.13 0.3 0.775 0.6]);
  set(hndl_win, 'currentaxes',h_axes);
  bar(eigenvalues);
  title('Eigenvalues');

  uiwait(hndl_win);
  firstEig = str2num(get(e_first, 'String'));
  lastEig = str2num(get(e_last, 'String'));

  % close the window
  close(hndl_win);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% See if the user has reduced the dimension enought

if lastEig > maxLastEig
  lastEig = maxLastEig;
  if b_verbose
    fprintf('Dimension reduced to %d due to the singularity of covariance matrix\n',...
           lastEig-firstEig+1);
  end
else
  % Reduce the dimensionality of the problem.
  if b_verbose
    if oldDimension == (lastEig - firstEig + 1)
      fprintf ('Dimension not reduced.\n');
    else
      fprintf ('Reducing dimension...\n');
    end
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drop the smaller eigenvalues
if lastEig < oldDimension
  lowerLimitValue = (eigenvalues(lastEig) + eigenvalues(lastEig + 1)) / 2;
else
  lowerLimitValue = eigenvalues(oldDimension) - 1;
end

lowerColumns = diag(D) > lowerLimitValue;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Drop the larger eigenvalues
if firstEig > 1
  higherLimitValue = (eigenvalues(firstEig - 1) + eigenvalues(firstEig)) / 2;
else
  higherLimitValue = eigenvalues(1) + 1;
end
higherColumns = diag(D) < higherLimitValue;

% Combine the results from above
selectedColumns = lowerColumns & higherColumns;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% print some info for the user
if b_verbose
  fprintf ('Selected [ %d ] dimensions.\n', sum (selectedColumns));
end
if sum (selectedColumns) ~= (lastEig - firstEig + 1),
  error ('Selected a wrong number of dimensions.');
end

if b_verbose
  fprintf ('Smallest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(lastEig));
  fprintf ('Largest remaining (non-zero) eigenvalue [ %g ]\n', eigenvalues(firstEig));
  fprintf ('Sum of removed eigenvalues [ %g ]\n', sum(diag(D) .* ...
    (~selectedColumns)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Select the colums which correspond to the desired range
% of eigenvalues.
E = selcol(E, selectedColumns);
D = selcol(selcol(D, selectedColumns)', selectedColumns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Some more information
if b_verbose
  sumAll=sum(eigenvalues);
  sumUsed=sum(diag(D));
  retained = (sumUsed / sumAll) * 100;
  fprintf('[ %g ] %% of (non-zero) eigenvalues retained.\n', retained);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newMatrix = selcol(oldMatrix, maskVector);

% newMatrix = selcol(oldMatrix, maskVector);
%
% Selects the columns of the matrix that marked by one in the given vector.
% The maskVector is a column vector.

% 15.3.1998

if size(maskVector, 1) ~= size(oldMatrix, 2),
  error ('The mask vector and matrix are of uncompatible size.');
end

numTaken = 0;

for i = 1 : size (maskVector, 1),
  if maskVector(i, 1) == 1,
    takingMask(1, numTaken + 1) = i;
    numTaken = numTaken + 1;
  end
end

newMatrix = oldMatrix(:, takingMask);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newVectors, meanValue] = remmean(vectors);
%REMMEAN - remove the mean from vectors
%
% [newVectors, meanValue] = remmean(vectors);
%
% Removes the mean of row vectors.
% Returns the new vectors and the mean.
%
% This function is needed by FASTICA and FASTICAG

% 24.8.1998
% Hugo Gävert

newVectors = zeros (size (vectors));
meanValue = mean (vectors')';
newVectors = vectors - meanValue * ones (1,size (vectors, 2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [newVectors, whiteningMatrix, dewhiteningMatrix] = whitenv ...
    (vectors, E, D, s_verbose);
%WHITENV - Whitenv vectors.
%
% [newVectors, whiteningMatrix, dewhiteningMatrix] = ...
%                               whitenv(vectors, E, D, verbose);
%
% Whitens the data (row vectors) and reduces dimension. Returns
% the whitened vectors (row vectors), whitening and dewhitening matrices.
%
% ARGUMENTS
%
% vectors       Data in row vectors.
% E             Eigenvector matrix from function 'pcamat'
% D             Diagonal eigenvalue matrix from function 'pcamat'
% verbose       Optional. Default is 'on'
%
% EXAMPLE
%       [E, D] = pcamat(vectors);
%       [nv, wm, dwm] = whitenv(vectors, E, D);
%
%
% This function is needed by FASTICA and FASTICAG
%
%   See also PCAMAT

% ========================================================
% Default value for 'verbose'
if nargin < 4, s_verbose = 'on'; end

% Check the optional parameter verbose;
switch lower(s_verbose)
 case 'on'
  b_verbose = 1;
 case 'off'
  b_verbose = 0;
 otherwise
  error(sprintf('Illegal value [ %s ] for parameter: ''verbose''\n', s_verbose));
end

% ========================================================
% Calculate the whitening and dewhitening matrices (these handle
% dimensionality simultaneously).
whiteningMatrix = inv (sqrt (D)) * E';
dewhiteningMatrix = E * sqrt (D);

% Project to the eigenvectors of the covariance matrix.
% Whiten the samples and reduce dimension simultaneously.
if b_verbose, fprintf ('Whitening...\n'); end
newVectors =  whiteningMatrix * vectors;

% ========================================================
% Just some security...
if ~isreal(newVectors)
  error ('Whitened vectors have imaginary values.');
end

% Print some information to user
if b_verbose
  fprintf ('Check: covariance differs from identity by [ %g ].\n', ...
    max (max (abs (cov (newVectors', 1) - eye (size (newVectors, 1))))));
end
