function [W,cost,isi] = icatb_iva_laplace_sos(X,varargin)
% [W,cost,shapePass,isi] = iva_mpe_decp_v2(X,varargin)
%
% Implementation of the independent vector analysis (IVA)
% algorithm using an correlated multivariate Laplace prior. This implementation also makes use of decoupling method to
% achieve 'Newton' like convergence. (see references below)
%
% Input:
% X - data observations from K data sets, i.e. X(:,:,k)=A(:,:,k)*S(:,:,k), 
% where A(:,:,k) is an N x N unknowned invertible mixing matrix and 
% S(:,:,k) is an N x T  matrix with the nth row corresponding to T samples 
% of the nth source in the kth dataset.  For IVA it is assumed that the 
% source is statistically independent of all the sources within a dataset 
% and exactly dependent on at most one source in each of the other 
% datasets. The data, X, is a 3-dimensional matrix of dimensions N x T x K.
% The latter enforces the assumption of an equal number of samples in each 
% dataset.
% 
% varargin - an optional input which is a structure containing the
% following elements. Note that default values are listed next to the 
% structure elements.
%
%    'whiten',true, ... % whitening is optional
%    'verbose',false, ... % verbose true enables print statements
%    'A',[], ... % true mixing matrices A, automatically sets verbose
%    'W_init',[], ... % initial estimates for demixing matrices in W
%    'maxIter',2*512, ... % max number of iterations
%    'terminationCriterion','ChangeInW', ... % criterion for terminating iterations: (ChangeInCost), ChangeInW
%    'termThreshold',1e-4, ... % termination threshold
%    'alpha0',0.1 ... % initial step size scaling
%
% Output:
% W - the estimated demixing matrices so that ideally 
% W(:,:,k)*A(:,:,k) = P*D(:,:,k) where P is any arbitrary permutation 
% matrix and D(:,:,k) is any diagonal invertible (scaling) matrix.  Note P 
% is common to all datasets; this is to indicate that the local permuation 
% ambiguity between dependent sources across datasets should ideally be 
% resolved by IVA.
%
% cost - the cost for each iteration
% 
% isi - joint inter-symbol-interference is available if user supplies true
% mixing matrices for computing a performance metric. Returns NaN if mixing
% matrix is not given as input.
%
%
% Example call:
% W=iva_mpe_decp_v2(X)
%
% To test code call function with no input arguments: iva_mpe_decp_v2
%
% Coded by Matthew Anderson (matt.anderson at umbc.edu)
% Modified by Suchita Bhinge (suchita1 at umbc.edu) on September 14 2018:
% Modification is in the score function in order to account for effect on
% large number of datasets
%
% References:
%
% [1] S. Bhinge, R. Mowakeaa, V. D. Calhoun, & T. Adali, "Extraction of time-varying spatio-temporal networks using parameter-tuned constrained IVA," Transaction on Medical Imaging, accepted for publication.
% [2] M. Anderson, T. Adali, & X.-L. Li,  "Joint Blind Source Separation of Multivariate Gaussian Sources: Algorithms and Performance Analysis," IEEE Trans. Signal Process., 2012, 60, 1672-1683
% [3] M. Anderson, G.-S. Fu, R. Phlypo, and T. Adali, "Independent vector analysis: Identification conditions and performance bounds," IEEE Trans. Signal Processing, vol. 62, pp. 4399-4410, Sep. 2014.

% Version 01 - 20120919 - Initial publication.
% Version 02 - 20140806 - Improved notation in comments.

%% Call internal test function when no arguments supplied
if nargin==0
   help iva_mpe_decp_v2
   test_iva_mpe_decp_v2
   return
end

%% Gather Options

% build default options structure
Params=struct( ...
   'whiten',true, ... % whitening is optional
   'gradProjection',false, ... % project gradient onto orthogonal direction
   'verbose',true, ... % verbose true enables print statements
   'A',[], ... % true mixing matrices A, automatically sets verbose
   'initW',[], ... % initial estimates for demixing matrices in W
   'maxIter',2*512, ... % max number of iterations
   'terminationCriterion','ChangeInW', ... % criterion for terminating iterations: ChangeInCost, (ChangeInW)
   'termThreshold',1e-4, ... % termination threshold
   'alpha0',1.0 ... % initial step size scaling
   );

% load in user supplied options
Params=getopt(Params,varargin{:});
alphaMin=Params.termThreshold; % alpha0 will max(alphaMin,alphaScale*alpha0)
alphaScale=0.9;
supplyA=~isempty(Params.A); % set to true if user has supplied true mixing matrices
outputISI=false;

%% Create cell versions
% Faster than using the more natural 3-dimensional version
if iscell(Params.A)
   K=length(Params.A);
   A=zeros([size(Params.A{1}) K]);
   for k=1:K
      A(:,:,k)=Params.A{k};
   end
   Params.A=A;
   clear A
end

[N,T,K]=size(X); % the input format insures all datasets have equal number of samples

if Params.whiten
   [X,V]=whiten(X);
else
    %Whitening matrix is identity matrix, needed for seed
    for k = 1 : K
        V(:,:,k) = diag(ones(1,N));
    end
end

Rx=cell(K,K);
for k1=1:K
   Rx{k1,k1}=1/T*(X(:,:,k1)*X(:,:,k1)');
   for k2=(k1+1):K
      Rx{k1,k2}=1/T*(X(:,:,k1)*X(:,:,k2)');
      Rx{k2,k1}=Rx{k1,k2}';
   end
end


%% Initialize W
if ~isempty(Params.initW)
   W=Params.initW;
   if size(W,3)==1 && size(W,1)==N && size(W,2)==N
      W=repmat(W,[1,1,K]);
   end
   if Params.whiten
      for k=1:K
         W(:,:,k)=W(:,:,k)/V(:,:,k);
      end
   end
else
   W=randn(N,N,K);
end

for k=1:K
   W(:,:,k)=vecnorm(W(:,:,k)')';
end
%% When mixing matrix A supplied
% verbose is set to true
% outputISI can be computed if requested
% A matrix is conditioned by V if data is whitened
if supplyA
   % only reason to supply A matrices is to display running performance
   %Params.verbose=true;
   if nargout>2
      outputISI=true;
   end
   if Params.whiten
      for k = 1:K
         Params.A(:,:,k) = V(:,:,k)*Params.A(:,:,k);
      end
   end
end

%% Initialize some local variables
cost=nan(1,Params.maxIter);
gammaShape=(K+1).^0.5; %.. Modified by Suchita Bhinge... September 14 2018
isi=nan(1,Params.maxIter);
Y=X*0;

%% Main Iteration Loop
for iter = 1:Params.maxIter
   termCriterion=0;
   
   % Current estimated sources
   for k=1:K
      Y(:,:,k)=W(:,:,k)*X(:,:,k);
   end
   
   % Some additional computations of performance via ISI when true A is supplied
   if supplyA
      [amari_avg_isi,amari_joint_isi]=bss_isi(W,Params.A);
      if outputISI
         isi(iter)=amari_joint_isi;
      end
   end
   
   W_old=W; % save current W as W_old
   [cost(iter),~,ft(iter),st(iter)]=comp_l_sos_cost(W,X);
   Q=0; R=0;
   
   for n=1:N      
      
      [hnk,Q,R]=decouple_trick(W,n,Q,R);
      
      yn=squeeze(Y(n,:,:))';
      % Efficient version of Ryn=yn*yn'/T;
      Ryn=eye(K); %
      for k1=1:K
         for k2=k1:K
            Ryn(k1,k2)=W(n,:,k1)*Rx{k1,k2}*W(n,:,k2)';
            if k1~=k2
               Ryn(k2,k1)=Ryn(k1,k2)';
            end
         end % k2
      end %k1
      
      %% Loop over each dataset
      for k=1:K
         % Derivative of cost function with respect to wnk
         invRyn=inv(Ryn);
         gipyn=dot(yn,invRyn*yn); %#ok<MINV> % since inverse already inverted, this is faster than Cov_n\yn
         
         phi=(invRyn(k,:)*yn).*gipyn.^(-0.5)*gammaShape; %.. Modified by Suchita Bhinge... September 14 2018
         dW= (X(:,:,k)*phi')/T - hnk(:,k)/(W(n,:,k)*hnk(:,k));
         
         if Params.gradProjection
            dW=vecnorm(dW - W(n,:,k)*dW*W(n,:,k)'); % non-colinear direction normalized
         end
         
         W(n,:,k)=vecnorm(W(n,:,k)' - Params.alpha0*dW)';
         yn(k,:)=W(n,:,k)*X(:,:,k);
         
         for kk=1:K
            Ryn(k,kk)=W(n,:,k)*Rx{k,kk}*W(n,:,kk)'; % = yn*yn'/T;
         end
         Ryn(:,k)=Ryn(k,:)';
      end % k
   end % n
   
   %% Calculate termination criterion
   switch lower(Params.terminationCriterion)
      case lower('ChangeInW')
         for k=1:K
            termCriterion = max(termCriterion,max(1-abs(diag(W_old(:,:,k)*W(:,:,k)'))));
         end % k
      case lower('ChangeInCost')
         if iter==1
            termCriterion=1;
         else
            termCriterion=abs(cost(iter-1)-cost(iter))/abs(cost(iter));
         end
      otherwise
         error('Unknown termination method.')
   end
   
   %% Check the termination condition
   if termCriterion < Params.termThreshold || iter == Params.maxIter
       
      break;
   elseif isnan(cost(iter))
      for k = 1:K
         W(:,:,k) = eye(N) + 0.1*randn(N);
      end
      if Params.verbose
         fprintf('\n W blowup, restart with new initial value.');
      end
   elseif iter>1 && cost(iter)>cost(iter-1)
      % see if this improves convergence
      Params.alpha0=max(alphaMin,alphaScale*Params.alpha0);
   end
   
   %% Display Iteration Information
   if Params.verbose
      if supplyA
         fprintf('\n Step %d: W change: %f, Cost: %f, Avg ISI: %f, Joint ISI: %f',  ...
            iter, termCriterion,cost(iter),amari_avg_isi,amari_joint_isi);
      else
         fprintf('\n Step %d: W change: %f, Cost: %f', iter, termCriterion,cost(iter));
      end
   end % options.verbose
end % iter

%% Finish Display
if iter==1 && Params.verbose
   if supplyA
      fprintf('\n Step %d: W change: %f, Cost: %f, Avg ISI: %f, Joint ISI: %f',  ...
         iter, termCriterion,cost(iter),amari_avg_isi,amari_joint_isi);
   else
      fprintf('\n Step %d: W change: %f, Cost: %f', iter, termCriterion,cost(iter));
   end
end % options.verbose

if Params.verbose
   fprintf('\n');
end

%% Clean-up Outputs
if Params.whiten
   for k=1:K
      W(:,:,k)=W(:,:,k)*V(:,:,k);
   end
end


%% Return

return;
% END OF iva_mpe_decp_v2

function test_iva_mpe_decp_v2
%% built-in test function for iva_mpe_decp_v2
disp('Executing internal test function for iva_mpe_decp_v2.')

%% Simulate input data
N=3;
K=50;
T=1e4;

S=zeros(N,T,K);
for n=1:N
   Z=randmv_laplace(K,T);
   S(n,:,:)=Z.';
end
A=randn(N,N,K);
X=S;
for k=1:K
   A(:,:,k)=vecnorm(A(:,:,k))';
   X(:,:,k)=A(:,:,k)*S(:,:,k);
end

%% Example call to iva_mpe_decp_v2
% with optional input, A
W=iva_mpe_decp_v2(X,'A',A);

%% Compute & display data
[~,isi]=bss_isi(W,A,S);
disp(['IVA-L-Decp ISI: ' num2str(isi)])
disp('IVA-L-Decp can "get stuck" in local minima, occurences of local minima seems to decrease with number of datasets.')


%% Visualize results
T1 = zeros(N,N);
for k = 1 : K
   Tk = W(:,:,k)*A(:,:,k);
   Tk = abs(Tk);
   for n = 1 : N
      Tk(n,:) = Tk(n,:)/max(abs(Tk(n,:)));
   end
   T1 = T1 + Tk/K;
end
P=zeros(N);
[~,imax]=max(T1);
ind=sub2ind([N N],1:N,imax);
P(ind)=1;
T1=P*T1;
figure
imagesc(1:N,1:N,T1)
caxis([0 1])

colormap('bone')
title('joint global matrix')
colorbar
disp('Ideally image is identity matrix.')
%% Return
return

function properties = getopt(properties,varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
%   No validation of the values is performed.
%%
% EXAMPLE:
%   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
%   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% would return:
%   properties =
%         zoom: 1
%       aspect: 0.7600
%        gamma: 1
%         file: 'mydata.dat'
%           bg: []
%
% Typical usage in a function:
%   properties = getopt(properties,varargin{:})

% Function from
% http://mathforum.org/epigone/comp.soft-sys.matlab/sloasmirsmon/bp0ndp$crq5@cui1.lmms.lmco.com

% dgleich
% 2003-11-19
% Added ability to pass a cell array of properties

if ~isempty(varargin) && (iscell(varargin{1}))
   varargin = varargin{1};
end;

% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii=1:length(varargin)
   arg = varargin{ii};
   if isempty(TargetField)
      if ~ischar(arg)
         error('Property names must be character strings');
      end
      %f = find(strcmp(prop_names, arg));
      if isempty(find(strcmp(prop_names, arg),1)) %length(f) == 0
         error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
      end
      TargetField = arg;
   else
      properties.(TargetField) = arg;
      TargetField = '';
   end
end
if ~isempty(TargetField)
   error('Property names and values must be specified in pairs.');
end


function mag=vecmag(vec,varargin)
% mag=vecmag(vec)
% or
% mag=vecmag(v1,v2,...,vN)
%
% Computes the vector 2-norm or magnitude of vector. vec has size n by m
% represents m vectors of length n (i.e. m column-vectors). Routine avoids
% potential mis-use of norm built-in function. Routine is faster than
% calling sqrt(dot(vec,vec)) -- but equivalent.
if nargin==1
   mag=sqrt(sum(vec.*conj(vec)));
else
   mag=vec.*conj(vec);
   for ii=1:length(varargin)
      mag=mag+varargin{ii}.*conj(varargin{ii});
   end
   mag=sqrt(mag);
end
return

function [uvec,mag]=vecnorm(vec)
% [vec,mag]=vecnorm(vec)
% Returns the vector normalized by 2-norm or magnitude of vector.
% vec has size n by m represents m vectors of length n (i.e. m
% column-vectors).
[n,m]=size(vec);
if n==1
   disp('vecnorm operates on column vectors, input appears to have dimension of 1')
end

uvec=zeros(n,m);
mag=vecmag(vec); % returns a 1 x m row vector
for ii=1:size(vec,1)
   uvec(ii,:)=vec(ii,:)./mag;
end
% Equivalent to: uvec=vec./repmat(mag,size(vec,1),1);

% Which implementation is optimal depends on optimality criterion (memory
% vs. speed), this version uses the former criterion.
return

function Y=randmv_laplace(d,T,varargin)
% Y=randmv_laplace(d,T)
% or
% Y=randmv_laplace(d,T,'lambda',lambd,'mu',mu,'Gamma',Gamm)
% where all additional parameters are optional with the following default
% values:
%
%  'lambda', 1, ... % exponential rate parameter, > 0
%  'mu', zeros(d,1), ... % mean vector
%  'Gamma', orth(rand(d)) ... % internal covariance structure, note det(Gamma)=1
%
% Generate T iid samples of the d-dimensional multivariate Laplace (ML)
% distribution, as given in "On the Multivariate Laplace Distribution", by
% Eltoft in IEEE Sig. Proc. Letters 2006.
%
% Note that a method for transforming an uncorrelated ML, Y, into a
% correlated ML, V, is given by Eltoft reference using V=A*Y+b, where A is d
% x d real-valued matrix and b is d x 1 real-valued vector, then V is ML
% with parameters specified in equations (14)-(16):
%
% lambda_new = lambda*abs(det(A))^(1/d)
% mu_new = A*mu + b
% Gamma_new = A*Gamma*A'*abs(det(A))^(-2/d)
%
% The special case of mu=b=0, Gamma=eye(d), and det(A)=1 is nice since
% then, lambda_new=lambda, mu_new=0, and Gamma_new=A*Gamma*A'.
%
% Coded by Matthew Anderson 9/15/2011.

global EXPRND_AVAILABLE RANDRAW_AVAILABLE

if isempty(EXPRND_AVAILABLE)
   EXPRND_AVAILABLE=true;
   if ~exist('exprnd.m','file')
      EXPRND_AVAILABLE=false;
   end
end
if ~EXPRND_AVAILABLE && isempty(RANDRAW_AVAILABLE)
   if exist('randraw.m','file')
      RANDRAW_AVAILABLE=true;
   else
      error('The exprnd function is unavailable, download randraw function via Matlab file exchange to proceed without the statistics toolbox.')
   end
end

if nargin==0
   test_randmv_laplace
   return
end

Params=struct( ...
   'lambda', 1, ... % exponential rate parameter, > 0
   'mu', zeros(d,1), ... % mean vector
   'Gamma', eye(d) ... % internal covariance structure, note det(Gamma)=1
   );
Params=getopt(Params,varargin{:});
if Params.lambda<0 || imag(Params.lambda)~=0
   error('Rate should be real-valued and greater than 0.')
end
Params.mu=Params.mu(:);
if length(Params.mu)~=d || any(imag(Params.lambda)~=0)
   error('Mean vector should be real-valued and correct dimension (d x 1).')
end
if size(Params.Gamma,1)~=d || size(Params.Gamma,2)~=d || ...
      abs(det(Params.Gamma)-1)>0.0001 || any(imag(Params.Gamma(:))~=0)
   error('Internal covariance structure needs to be real-valued square matrix with a determinant of one.')
end

X = randn(d,T);
if EXPRND_AVAILABLE
   Z = sqrt(exprnd(1/Params.lambda,[1 T]));
else %if RANDRAW_AVAILABLE
   Z = sqrt(randraw('exp', Params.lambda, [1 T]));
end
Y = bsxfun(@plus,Params.mu,bsxfun(@times,Z,sqrtm(Params.Gamma)*X)); % equation (6)

return

function test_randmv_laplace
help randmv_laplace
disp('Executing internal test function for randmv_laplace.')
T=1e6;
dim=2;
Nbins=[100 100];

Gamm=[1 0.5; 0.5 1];
Gamm=Gamm/sqrt(det(Gamm));
Gamm=eye(dim);
lambd=1;
mu=[0 2]';
Z=randmv_laplace(dim,T,'Gamma',Gamm,'lambda',lambd,'mu',mu);

figure
[hist_out,hist_bins]=hist3(Z','Nbins',Nbins);
imagesc(hist_bins{1},hist_bins{2},hist_out')
axis xy; colorbar
xlabel('Z(1)')
ylabel('Z(2)')
title(['Histogram of Multivariate Laplace'])

[z1,z2]=meshgrid(hist_bins{1},hist_bins{2});
f=pdf_mvlaplace([z1(:)'; z2(:)'],lambd,Gamm,mu);

figure
f=reshape(f,length(hist_bins{1}),length(hist_bins{2}));
imagesc(hist_bins{1},hist_bins{2}, f)
axis xy; colorbar
xlabel('Z(1)')
ylabel('Z(2)')
title(['pdf of Multivariate Laplace'])

return

function [z,V,U]=whiten(x)
% [z,V,U]=whiten(x)
%
% Whitens the data vector so that E{zz'}=I, where z=V*x.

if ~iscell(x)
   [N,T,K]=size(x);
   if K==1
      % Step 1. Center the data.
      x=bsxfun(@minus,x,mean(x,2));
      
      % Step 2. Form MLE of data covariance.
      covar=x*x'/T;
      
      % Step 3. Eigen decomposition of covariance.
      [eigvec, eigval] = eig (covar);
      
      % Step 4. Forming whitening transformation.
      V=sqrt(eigval) \ eigvec';
      U=eigvec * sqrt(eigval);
      
      % Step 5. Form whitened data
      z=V*x;
   else
      K=size(x,3);
      z=zeros(N,T,K);
      V=zeros(N,N,K);
      U=zeros(N,N,K);
      for k=1:K
         % Step 1. Center the data.
         xk=bsxfun(@minus,x(:,:,k),mean(x(:,:,k),2));
         
         % Step 2. Form MLE of data covariance.
         covar=xk*xk'/T;
         
         % Step 3. Eigen decomposition of covariance.
         [eigvec, eigval] = eig (covar);
         
         % Step 4. Forming whitening transformation.
         V(:,:,k)=sqrt(eigval) \ eigvec';
         U(:,:,k)=eigvec * sqrt(eigval);
         
         % Step 5. Form whitened data
         z(:,:,k)=V(:,:,k)*xk;
      end % k
   end % K>1
else % x is cell
   K=numel(x);
   sizex=size(x);
   V=cell(sizex);
   U=cell(sizex);
   z=cell(sizex);
   for k=1:K
      T=size(x{k},2);
      % Step 1. Center the data.
      xk=bsxfun(@minus,x{k},mean(x{k},2));
      
      % Step 2. Form MLE of data covariance.
      covar=xk*xk'/T;
      
      % Step 3. Eigen decomposition of covariance.
      [eigvec, eigval] = eig (covar);
      
      % Step 4. Forming whitening transformation.
      V{k}=sqrt(eigval) \ eigvec';
      U{k}=eigvec * sqrt(eigval);
      
      % Step 5. Form whitened data
      z{k}=V{k}*xk;
   end % k
end

%%
return

function [isi,isiGrp,success,G]=bss_isi(W,A,s,Nuse)
% Non-cell inputs:
% isi=bss_isi(W,A) - user provides W & A where x=A*s, y=W*x=W*A*s
% isi=bss_isi(W,A,s) - user provides W, A, & s
%
% Cell array of matrices:
% [isi,isiGrp]=bss_isi(W,A) - W & A are cell array of matrices
% [isi,isiGrp]=bss_isi(W,A,s) - W, A, & s are cell arrays
%
% 3-d Matrices:
% [isi,isiGrp]=bss_isi(W,A) - W is NxMxK and A is MxNxK
% [isi,isiGrp]=bss_isi(W,A,s) - S is NxTxK (N=#sources, M=#sensors, K=#datasets)
%
% Measure of quality of separation for blind source separation algorithms.
% W is the estimated demixing matrix and A is the true mixing matrix.  It should be noted
% that rows of the mixing matrix should be scaled by the necessary constants to have each
% source have unity variance and accordingly each row of the demixing matrix should be
% scaled such that each estimated source has unity variance.
%
% ISI is the performance index given in Complex-valued ICA using second order statisitcs
% Proceedings of the 2004 14th IEEE Signal Processing Society Workshop, 2004, 183-192
%
% Normalized performance index (Amari Index) is given in Choi, S.; Cichocki, A.; Zhang, L.
% & Amari, S. Approximate maximum likelihood source separation using the natural gradient
% Wireless Communications, 2001. (SPAWC '01). 2001 IEEE Third Workshop on Signal
% Processing Advances in, 2001, 235-238.
%
% Note that A is p x M, where p is the number of sensors and M is the number of signals
% and W is N x p, where N is the number of estimated signals.  Ideally M=N but this is not
% guaranteed.  So if N > M, the algorithm has estimated more sources than it "should", and
% if M < N the algorithm has not found all of the sources.  This meaning of this metric is
% not well defined when averaging over cases where N is changing from trial to trial or
% algorithm to algorithm.

% Some examples to consider
% isi=bss_isi(eye(n),eye(n))=0
%
% isi=bss_isi([1 0 0; 0 1 0],eye(3))=NaN
%


% Should ideally be a permutation matrix with only one non-zero entry in any row or
% column so that isi=0 is optimal.

% generalized permutation invariant flag (default=false), only used when nargin<3
gen_perm_inv_flag=false;
success=true;

Wcell=iscell(W);
if nargin<2
   Acell=false;
else
   Acell=iscell(A);
end
if ~Wcell && ~Acell
   if ndims(W)==2 && ndims(A)==2
      if nargin==2
         % isi=bss_isi(W,A) - user provides W & A
         
         % Traditional Metric, user provided W & A separately
         G=W*A;
         [N,M]=size(G);
         Gabs=abs(G);
         if gen_perm_inv_flag
            % normalization by row
            max_G=max(Gabs,[],2);
            Gabs=repmat(1./max_G,1,size(G,2)).*Gabs;
         end
      elseif nargin==3
         % Equalize energy associated with each estimated source and true
         % source.
         %
         % y=W*A*s;
         % snorm=D*s; where snorm has unit variance: D=diag(1./std(s,0,2))
         % Thus: y=W*A*inv(D)*snorm
         % ynorm=U*y; where ynorm has unit variance: U=diag(1./std(y,0,2))
         % Thus: ynorm=U*W*A*inv(D)*snorm=G*snorm and G=U*W*A*inv(D)
         
         y=W*A*s;
         D=diag(1./std(s,0,2));
         U=diag(1./std(y,0,2));
         G=U*W*A/D; % A*inv(D)
         [N,M]=size(G);
         Gabs=abs(G);
      else
         error('Not acceptable.')
      end
      
      isi=0;
      for n=1:N
         isi=isi+sum(Gabs(n,:))/max(Gabs(n,:))-1;
      end
      for m=1:M
         isi=isi+sum(Gabs(:,m))/max(Gabs(:,m))-1;
      end
      isi=isi/(2*N*(N-1));
      isiGrp=NaN;
      success=NaN;
   elseif ndims(W)==3 && ndims(A)==3
      % IVA/GroupICA/MCCA Metrics
      % For this we want to average over the K groups as well as provide the additional
      % measure of solution to local permutation ambiguity (achieved by averaging the K
      % demixing-mixing matrices and then computing the ISI of this matrix).
      [N,M,K]=size(W);
      if M~=N
         error('This more general case has not been considered here.')
      end
      L=M;
      
      isi=0;
      GabsTotal=zeros(N,M);
      G=zeros(N,M,K);
      for k=1:K
         if nargin<=2
            % Traditional Metric, user provided W & A separately
            Gk=W(:,:,k)*A(:,:,k);
            Gabs=abs(Gk);
            if gen_perm_inv_flag
               % normalization by row
               max_G=max(Gabs,[],2);
               Gabs=repmat(1./max_G,1,size(Gabs,2)).*Gabs;
            end
         else %if nargin==3
            % Equalize energy associated with each estimated source and true
            % source.
            %
            % y=W*A*s;
            % snorm=D*s; where snorm has unit variance: D=diag(1./std(s,0,2))
            % Thus: y=W*A*inv(D)*snorm
            % ynorm=U*y; where ynorm has unit variance: U=diag(1./std(y,0,2))
            % Thus: ynorm=U*W*A*inv(D)*snorm=G*snorm and G=U*W*A*inv(D)
            yk=W(:,:,k)*A(:,:,k)*s(:,:,k);
            Dk=diag(1./std(s(:,:,k),0,2));
            Uk=diag(1./std(yk,0,2));
            Gk=Uk*W(:,:,k)*A(:,:,k)/Dk;
            
            Gabs=abs(Gk);
         end
         G(:,:,k)=Gk;
         
         if nargin>=4
            Np=Nuse;
            Mp=Nuse;
            Lp=Nuse;
         else
            Np=N;
            Mp=M;
            Lp=L;
         end
         
         % determine if G is success by making sure that the location of maximum magnitude in
         % each row is unique.
         if k==1
            [~,colMaxG]=max(Gabs,[],2);
            if length(unique(colMaxG))~=Np
               % solution is failure in strictest sense
               success=false;
            end
         else
            [~,colMaxG_k]=max(Gabs,[],2);
            if ~all(colMaxG_k==colMaxG)
               % solution is failure in strictest sense
               success=false;
            end
         end
         
         GabsTotal=GabsTotal+Gabs;
         
         for n=1:Np
            isi=isi+sum(Gabs(n,:))/max(Gabs(n,:))-1;
         end
         for m=1:Mp
            isi=isi+sum(Gabs(:,m))/max(Gabs(:,m))-1;
         end
      end
      isi=isi/(2*Np*(Np-1)*K);
      
      Gabs=GabsTotal;
      if gen_perm_inv_flag
         % normalization by row
         max_G=max(Gabs,[],2);
         Gabs=repmat(1./max_G,1,size(Gabs,2)).*Gabs;
      end
      %       figure; imagesc(Gabs); colormap('bone'); colorbar
      isiGrp=0;
      for n=1:Np
         isiGrp=isiGrp+sum(Gabs(n,:))/max(Gabs(n,:))-1;
      end
      for m=1:Mp
         isiGrp=isiGrp+sum(Gabs(:,m))/max(Gabs(:,m))-1;
      end
      isiGrp=isiGrp/(2*Lp*(Lp-1));
   else
      error('Need inputs to all be of either dimension 2 or 3')
   end
elseif Wcell && Acell
   % IVA/GroupICA/MCCA Metrics
   % For this we want to average over the K groups as well as provide the additional
   % measure of solution to local permutation ambiguity (achieved by averaging the K
   % demixing-mixing matrices and then computing the ISI of this matrix).
   
   K=length(W);
   N=0; M=0;
   Nlist=zeros(K,1);
   for k=1:K
      Nlist(k)=size(W{k},1);
      N=max(size(W{k},1),N);
      M=max(size(A{k},2),M);
   end
   commonSources=false; % limits the ISI to first min(Nlist) sources
   if M~=N
      error('This more general case has not been considered here.')
   end
   L=M;
   
   % To make life easier below lets sort the datasets to have largest
   % dataset be in k=1 and smallest at k=K;
   [Nlist,isort]=sort(Nlist,'descend');
   W=W(isort);
   A=A(isort);
   if nargin > 2
      s=s(isort);
   end
   G=cell(K,1);
   isi=0;
   if commonSources
      minN=min(Nlist);
      GabsTotal=zeros(minN);
      Gcount=zeros(minN);
   else
      GabsTotal=zeros(N,M);
      Gcount=zeros(N,M);
   end
   for k=1:K
      if nargin==2
         % Traditional Metric, user provided W & A separately
         G{k}=W{k}*A{k};
         Gabs=abs(G{k});
         if gen_perm_inv_flag
            % normalization by row
            max_G=max(Gabs,[],2);
            Gabs=repmat(1./max_G,1,size(Gabs,2)).*Gabs;
         end
      elseif nargin>=3
         % Equalize energy associated with each estimated source and true
         % source.
         %
         % y=W*A*s;
         % snorm=D*s; where snorm has unit variance: D=diag(1./std(s,0,2))
         % Thus: y=W*A*inv(D)*snorm
         % ynorm=U*y; where ynorm has unit variance: U=diag(1./std(y,0,2))
         % Thus: ynorm=U*W*A*inv(D)*snorm=G*snorm and G=U*W*A*inv(D)
         yk=W{k}*A{k}*s{k};
         Dk=diag(1./std(s{k},0,2));
         Uk=diag(1./std(yk,0,2));
         G{k}=Uk*W{k}*A{k}/Dk;
         
         Gabs=abs(G{k});
      else
         error('Not acceptable.')
      end
      
      if commonSources
         Nk=minN;
         Gabs=Gabs(1:Nk,1:Nk);
      elseif nargin>=4
         commonSources=true;
         Nk=Nuse;
         minN=Nk;
      else
         Nk=Nlist(k);
      end
      
      if k==1
         [~,colMaxG]=max(Gabs(1:Nk,1:Nk),[],2);
         if length(unique(colMaxG))~=Nk
            % solution is a failure in a strict sense
            success=false;
         end
      elseif success
         if nargin>=4
            [~,colMaxG_k]=max(Gabs(1:Nk,1:Nk),[],2);
         else
            [~,colMaxG_k]=max(Gabs,[],2);
         end
         if ~all(colMaxG_k==colMaxG(1:Nk))
            % solution is a failure in a strict sense
            success=false;
         end
      end
      
      if nargin>=4
         GabsTotal(1:Nk,1:Nk)=GabsTotal(1:Nk,1:Nk)+Gabs(1:Nk,1:Nk);
      else
         GabsTotal(1:Nk,1:Nk)=GabsTotal(1:Nk,1:Nk)+Gabs;
      end
      Gcount(1:Nk,1:Nk)=Gcount(1:Nk,1:Nk)+1;
      for n=1:Nk
         isi=isi+sum(Gabs(n,:))/max(Gabs(n,:))-1;
      end
      for m=1:Nk
         isi=isi+sum(Gabs(:,m))/max(Gabs(:,m))-1;
      end
      isi=isi/(2*Nk*(Nk-1));
   end
   
   if commonSources
      Gabs=GabsTotal;
   else
      Gabs=GabsTotal./Gcount;
   end
   % normalize entries into Gabs by the number of datasets
   % contribute to each entry
   
   if gen_perm_inv_flag
      % normalization by row
      max_G=max(Gabs,[],2);
      Gabs=repmat(1./max_G,1,size(Gabs,2)).*Gabs;
   end
   isiGrp=0;
   
   if commonSources
      for n=1:minN
         isiGrp=isiGrp+sum(Gabs(n,1:minN))/max(Gabs(n,1:minN))-1;
      end
      for m=1:minN
         isiGrp=isiGrp+sum(Gabs(1:minN,m))/max(Gabs(1:minN,m))-1;
      end
      isiGrp=isiGrp/(2*minN*(minN-1));
   else
      for n=1:Nk
         isiGrp=isiGrp+sum(Gabs(n,:))/max(Gabs(n,:))-1;
      end
      for m=1:Nk
         isiGrp=isiGrp+sum(Gabs(:,m))/max(Gabs(:,m))-1;
      end
      isiGrp=isiGrp/(2*L*(L-1));
   end
   
else
   % Have not handled when W is cell and A is single matrix or vice-versa.  Former makes
   % sense when you want performance of multiple algorithms for one mixing matrix, while
   % purpose of latter is unclear.
end

return

function [h,invQ,R]=decouple_trick(W,n,invQ,R)
% h=decouple_trick(W,n)
% h=decouple_trick(W,n,1)
% [h,invQ]=decouple_trick(W,n,invQ)
% [h,Q,R]=decouple_trick(W,n,Q,R)
%
% Computes the h vector for the decoupling trick [1] of the nth row of W. W
% can be K 'stacked' square matrices, i.e., W has dimensions N x N x K.
% The output vector h will be formatted as an N x K matrix.  There are many
% possible methods for computing h.  This routine provides four different
% (but of course related) methods depending on the arguments used.
%
% Method 1:
% h=decouple_trick(W,n)
% h=decouple_trick(W,n,0)
% -Both calls above will result in the same algorithm, namely the QR
% algorithm is used to compute h.
%
% Method 2:
% h=decouple_trick(W,n,~), where ~ is anything
% -Calls the projection method.
%
% Method 3:
% [h,invQ]=decouple_trick(W,n,invQ)
% -If two output arguments are specified then the recursive algorithm
% described in [2].  It is assumed that the decoupling will be performed in
% sequence see the demo subfunction for details.
% An example call sequence:
%  [h1,invQ]=decouple_trick(W,1);
%  [h2,invQ]=decouple_trick(W,2,invQ);
%
% Method 4:
% [h,Q,R]=decouple_trick(W,n,Q,R)
% -If three output arguments are specified then a recursive QR algorithm is
% used to compute h.
% An example call sequence:
%  [h1,Q,R]=decouple_trick(W,1);
%  [h2,Q,R]=decouple_trick(W,2,Q,R);
%
% See the subfunction demo_decoupling_trick for more examples.  The demo
% can be executed by calling decouple_trick with no arguments, provides a
% way to compare the speed and determine the accuracy of all four
% approaches.
%
% Note that methods 2 & 3 do not normalize h to be a unit vector.  For
% optimization this is usually not of interest.  If it is then set the
% variable boolNormalize to true.
%
% Main References:
% [1] X.-L. Li & X.-D. Zhang, "Nonorthogonal Joint Diagonalization Free of Degenerate Solution," IEEE Trans. Signal Process., 2007, 55, 1803-1814
% [2] X.-L. Li & T. Adali, "Independent component analysis by entropy bound minimization," IEEE Trans. Signal Process., 2010, 58, 5151-5164
%
% Coded by Matthew Anderson (matt dot anderson at umbc dot edu)

% Version 01 - 20120919 - Initial publication


if nargin==0
   help decouple_trick
   demo_decouple_trick
   return
end
if nargin==1
   help decouple_trick
   error('Not enough inputs -- see displayed help above.')
end
[M,N,K]=size(W);
if M~=N
   error('Assuming W is square matrix.')
end
h=zeros(N,K);

% enables an additional computation that is usually not necessary if the
% derivative is of  interest, it is only necessary so that sqrt(det(W*W'))
% = sqrt(det(Wtilde*Wtilde'))*abs(w'*h) holds.  Furthermore, it is only
% necessary when using the recursive or projection methods.
%
% a user might wish to enable the calculation by setting the quantity below
% to true
boolNormalize=false;

if nargout==3
   % use QR recursive method
   % [h,Qnew,Rnew]=decouple_trick(W,n,Qold,Rold)
   if n==1
      invQ=zeros(N,N,K);
      R=zeros(N,N-1,K);
   end
   for k=1:K
      if n==1
         Wtilde=W(2:N,:,k);
         [invQ(:,:,k),R(:,:,k)]=qr(Wtilde');
      else
         n_last=n-1;
         e_last = zeros(N-1,1);
         e_last(n_last) = 1;
         [invQ(:,:,k),R(:,:,k)]=qrupdate(invQ(:,:,k),R(:,:,k),-W(n,:,k)',e_last);
         [invQ(:,:,k),R(:,:,k)]=qrupdate(invQ(:,:,k),R(:,:,k),W(n_last,:,k)',e_last);
      end
      h(:,k)=invQ(:,end,k); % h should be orthogonal to W(nout,:,k)'
   end
elseif nargout==2
   % use recursive method
   % [h,invQ]=decouple_trick(W,n,invQ), for any value of n=1, ..., N
   % [h,invQ]=decouple_trick(W,1), when n=1
   
   if n==1
      invQ=zeros(N-1,N-1,K);
   end
   % Implement a faster approach to calculating h.
   for k=1:K
      if n==1
         Wtilde=W(2:N,:,k);
         invQ(:,:,k)=inv(Wtilde*Wtilde');
      else
         if nargin<3
            help decouple_trick
            error('Need to supply invQ for recursive approach.')
         end
         [Mq,Nq,Kq]=size(invQ);
         if Mq~=(N-1) || Nq~=(N-1) || Kq~=K
            help decouple_trick
            error('Input invQ does not have the expected dimensions.')
         end
         n_last=n-1;
         Wtilde_last=W([(1:n_last-1) (n_last+1:N)],:,k);
         w_last=W(n_last,:,k)';
         w_current=W(n,:,k)';
         c = Wtilde_last*(w_last - w_current);
         c(n_last) = 0.5*( w_last'*w_last - w_current'*w_current );
         %e_last = zeros(N-1,1);
         %e_last(n_last) = 1;
         temp1 = invQ(:,:,k)*c;
         temp2 = invQ(:,n_last,k);
         inv_Q_plus = invQ(:,:,k) - temp1*temp2'/(1+temp1(n_last));
         
         temp1 = inv_Q_plus'*c;
         temp2 = inv_Q_plus(:,n_last);
         invQ(:,:,k) = inv_Q_plus - temp2*temp1'/(1+c'*temp2);
         % inv_Q is Hermitian
         invQ(:,:,k) = (invQ(:,:,k)+invQ(:,:,k)')/2;
      end
      
      temp1 = randn(N, 1);
      Wtilde = W([(1:n-1) (n+1:N)],:,k);
      h(:,k) = temp1 - Wtilde'*invQ(:,:,k)*Wtilde*temp1;
   end
   if boolNormalize
      h=vecnorm(h);
   end
elseif nargin==2 || invQ==0
   % use (default) QR approach
   % h=decouple_trick(W,n)
   % h=decouple_trick(W,n,0)
   for k=1:K
      [Q,~]=qr(W([(1:n-1) (n+1:N)],:,k)');
      h(:,k)=Q(:,end); % h should be orthogonal to W(nout,:,k)'
   end
else % use projection method
   % h=decouple_trick(W,n,~), ~ is anything
   for k=1:K
      temp1 = randn(N, 1);
      Wtilde = W([(1:n-1) (n+1:N)],:,k);
      h(:,k) = temp1 - Wtilde'*((Wtilde*Wtilde')\Wtilde)*temp1;
   end
   if boolNormalize
      h=vecnorm(h);
   end
end

return

function [cost,y,firstterm,secterm]=comp_l_sos_cost(W,X,const_log,uncorr)


if nargin<4
   uncorr=false;
   if nargin<3
      const_log=[];
   end
end
if iscell(W)
   error('not done for cells.')
   [N,T]=size(X{1});
   K=length(X);
   y=X;
   y_all=zeros(K,N,T);
   cost=0;

   for k = 1:K
      y{k} = W{k}*X{k};
      y_all(k,:,:) = y{k};  % subj X comp X T
      cost=cost-log(abs(det(W{k})));
   end

   for n=1:N
      y_comp=squeeze(y_all(:,n,:)); % subj X T
      gip=dot(y_comp,y_comp);
      dcost=mean(sqrt(gip));
      cost=cost+dcost;
   end
else
   scalesources=false;
   [N,T,K]=size(X);
   Y=X;
   cost=0;
   xi_shape=K+1;    %.. Modified by Suchita Bhinge... September 14 2018
   if isempty(const_log)
      const_Ck=0.5.*gamma(K/2).*xi_shape.^(K/2)./ ...
         (pi^(K/2)*gamma(K));      
      const_log=-log(const_Ck);
   end

   if scalesources
      for k = 1:K
         Y(:,:,k) = W(:,:,k)*X(:,:,k);
         ypower=diag(Y(:,:,k)*Y(:,:,k)')/T;
         W(:,:,k) = W(:,:,k)./repmat(sqrt(ypower),1,N);
      end
   end

   for k = 1:K
      Y(:,:,k) = W(:,:,k)*X(:,:,k);
      cost=cost-log(abs(det(W(:,:,k))));
   end % K
   
   firstterm = 0;
   secterm = 0;
   for n=1:N
      yn=shiftdim(Y(n,:,:)).'; % K x T
      if uncorr
         gip=dot(yn,yn);
         dcost=-K*log(rhat)+rhat*mean(sqrt(gip));
      else
         CovN=yn*yn'/T;
         gip=dot(yn,CovN\yn);
         dcost= (const_log + 0.5*log(det(CovN))) + xi_shape^0.5*mean(gip.^0.5);  %.. Modified by Suchita Bhinge... September 14 2018
      end
      secterm = secterm + (0.5*log(det(CovN)));
      firstterm = firstterm + xi_shape^0.5*mean(gip.^0.5);
      cost=cost+dcost;
   end
   y=Y;

   if nargout>5     
      error('Not done yet.')
      bool_numerical_derivative=false;
      cost_derivative=zeros(N,K,N); % derivative is computed assuming that each set of SCV demixing vectors are constant.
      if bool_numerical_derivative
         cost_derivative_numerical=zeros(N,K,N);
      end
      for n=1:N
         yn=shiftdim(Y(n,:,:)).'; % K x T
         SigN=yn*yn'/T;
         invSigN=inv(SigN);
         gip=dot(yn,invSigN*yn);         
         for k=1:K
            xx=bsxfun(@times,X(:,:,k),1./sqrt(gip));

            wnk=W(n,:,k); wnk=wnk(:);

            cost_derivative(n,k,:)=sqrt(K+1)*mean(bsxfun(@times,invSigN(k,:)*yn,xx),2);
            if 1==1
               % Find direction orthogonal to space spanned by excluded components demixing
               % matrix via QR decomposition
               [Q,R]=qr(W(nout,:,k)'); %#ok<NASGU>
               hnk=Q(:,end); % hn should be orthogonal to W(nout,:,k)'
               cost_derivative(n,k,:) = squeeze(cost_derivative(n,k,:)) - hnk/(hnk'*wnk);
            end

            if bool_numerical_derivative
               dd=0.00001;
               cost_d_nk=zeros(N,1);
               for nn=1:N
                  wnkplus=wnk;
                  wnkplus(nn)=wnkplus(nn)+dd;
                  wnkplus=vecnorm(wnkplus);
                  Wplus=W;
                  Wplus(n,:,k)=wnkplus';
                  costplus=comp_multivarlaplace_corr_cost(Wplus,X);

                  wnkminus=wnk;
                  wnkminus(nn)=wnkminus(nn)-dd;
                  wnkminus=vecnorm(wnkminus);
                  Wminus=W;
                  Wminus(n,:,k)=wnkminus';
                  costminus=comp_multivarlaplace_corr_cost(Wminus,X);

                  cost_d_nk(nn)=(costplus-costminus)/(2*dd);

               end % nn=1:N
               cost_derivative_numerical(n,k,:)=cost_d_nk;

               alph=0.001;
               dw=vecnorm(squeeze(cost_derivative(n,k,:)));
               dw_theory=vecnorm(dw-wnk'*dw*wnk);
               wnew_theory = vecnorm(wnk - alph*dw_theory);
               W_theory=W; W_theory(n,:,k)=wnew_theory;
               cost_theory=comp_multivarlaplace_corr_cost(W_theory,X);

               dw=vecnorm(squeeze(cost_derivative_numerical(n,k,:)));
               dw_numerical=vecnorm(dw-wnk'*dw*wnk);
               wnew_numerical = vecnorm(wnk - alph*dw_numerical);
               W_numerical=W; W_numerical(n,:,k)=wnew_numerical;
               cost_numerical=comp_multivarlaplace_corr_cost(W_numerical,X);

               disp(['n,k=' int2str(n) ', ' int2str(k) ': ' num2str(cost-cost_theory) ', ' ...
                  num2str(cost-cost_numerical)])
               disp([num2str(vecnorm(dw_numerical)'*vecnorm(dw_theory))])
            end % bool_numerical_derivative
         end % k=1:K
      end % n=1:N
      asdf=1;
   end % nargout>2
end
return