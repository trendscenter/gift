function [W,cost_function] = icatb_iva_laplace_sos(X,varargin)
% 11/10/25 Cyrus at TReNDs
% 11/10/25, Replaced function based upon update from Erdem at Adali lab
% 11/10/25, Previous icatb_iva_laplace_sos.m found in git records
% Trung's implementation - 02/2024
%
% Implementation of IVA with the multivariate laplacian distribution.
%
% For a general description of the algorithm and its relationship with others,
% see http://mlsp.umbc.edu/jointBSS_introduction.html
%
%%% Input:
% X - data observations from K data sets, i.e. X{k}=A{k}*S{k}, where A{k}
% is an N x N unknowned invertible mixing matrix and S{k} is N x T  matrix
% with the nth row corresponding to T samples of the nth source in the kth
% dataset.  For IVA it is assumed that the source is statistically
% independent of all the sources within a dataset and exactly dependent on
% at most one source in each of the other datasets. The data, X, is a 
% 3-dimensional matrix of dimensions N x T x K. 
% 
%%% Output:
% W - the estimated demixing matrices so that ideally W{k}*A{k} = P*D{k}
% where P is any arbitrary permutation matrix and D{k} is any diagonal
% invertible (scaling) matrix.  Note P is common to all datasets; this is
% to indicate that the local permuation ambiguity between dependent sources
% across datasets should ideally be resolved by IVA.
%
% cost_function - the objective function for each iteration
%
%%% Required functions in this package:
% decouple_trick.m
% getopt.m
% logdet.m
% pca_whitening.m
%
%%% Example call:
% W = iva_l_sos(X, 'initW', W0, 'minIter', 100)
% Optional input pairs and default values are given in Params
%
%%% References
% [1]- S. Bhinge, R. Mowakeaa, V.D. Calhoun, T. Adali, "Extraction of time-varying 
%    spatio-temporal networks using parameter-tuned constrained IVA." IEEE 
%    Transactions on Medical Imaging, 2019, vol. 38, no. 7, 1715-1725 
% [2] - Mowakeaa, R., Boukouvalas, Z., Long, Q., & Adali, T. (2019). "IVA using
%    complex multivariate GGD: application to fMRI analysis."
%    Multidimensional Systems and Signal Processing, 1-20.
% [3] - E. E. Kumbasar, H. Yang, V. Trung, V. D. Calhoun, & T. Adali, 
%    "Fusion of Multitask fMRI Data with Constrained Independent Vector 
%    Analysis,"2025 59th Annual Conference on Information Sciences and Systems 
%    (CISS), Baltimore, MD, USA, 2025, pp. 1-6.
% [3] - 
%%% Copyright (C) 2024 MLSP Lab
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% Please make sure you include the references in any related work.
% 
% A copy of the GNU General Public License can be found at
% <https://mlsp.umbc.edu/resources.html>.
% If not, see <https://www.gnu.org/licenses/>.


% Version 01 - 20120919 - Initial publication.
% Version 02 - 20140806 - Improved notation in comments.

%% Call internal test function when no arguments supplied
if nargin==0
   help iva_mpe_decp_v2
   test_iva_mpe_decp_v2
   return
end

    %% Gather Options
    [N, T, K] = size(X);
    
    %% Gather options
    
    % build default options structure
    Params = struct( ...
        'whiten', true, ... % whitening is optional
        'initW', [], ... % initial estimates for demixing matrices in W
        'maxIter', 1000, ... % max number of iterations
        'minIter', 300, ... % min number of iterations
        'termThreshold', 1e-6, ... % termination threshold based on change in W
        'alpha0', 1, ... % initial step size
        'stepsizeDecrease', 0.95, ... % decrement scale in step size if cost function increases
        'minStepsize', 1e-6, ... % minimum step size
        'checkpoints', false, ... % periodically save current results
        'updateCOV', true, ... % update the covariance estimate for every k
        'sortSCVs', true, ... % reoder SCVs from most to least ill-conditioned
        'checkpointID', '' ... % for naming the checkpoint files
    );
    
    % load in user supplied options
    Params = getopt(Params, varargin{:});
    beta = 0.5;
    stepsize = Params.alpha0;
    
    
    % compute Kotz constant while avoiding numerical error when K is large
    % and/or beta is small
    % COV = p_Kotz * Sigma
    log_p_Kotz = log(2) / beta + gammaln((K + 2)/2/beta) - log(K) - gammaln(K/2/beta);
    p_Kotz = exp(log_p_Kotz);
    
    if Params.whiten
        [X, V] = pca_whitening(X);
    end
    
    % Pre-compute correlation in X
    Rx = cell(K, K);
    for k = 1:K
        Rx{k, k} = 1 / T * (X(:, :, k) * X(:, :, k)');
        for k2 = (k + 1):K
            Rx{k, k2} = 1 / T * (X(:, :, k) * X(:, :, k2)');
            Rx{k2, k} = Rx{k, k2}';
        end
    end
    
    % Initialize W
    if ~isempty(Params.initW)
        W = Params.initW;
        if Params.whiten
            for k = 1:K
                W(:, :, k) = W(:, :, k) / V(:, :, k);
            end
        end
    else
        W = randn(N, N, K);
    end
    % make sure rows of W is normalized
    for k = 1:K
        W(:, :, k) = vecnorm(W(:, :, k));
    end
    
    % Initialize source-related quantities based on initial W
    COV_N_est = zeros(K, K, N);
    Y = zeros(size(X));
    for k = 1:K
        Y(:, :, k) = W(:, :, k) * X(:, :, k);
    
        for n = 1:N
            for k2 = k:K
                COV_N_est(k, k2, n) = W(n, :, k) * Rx{k, k2} * W(n, :, k2)';
                if k < k2
                    COV_N_est(k2, k, n) = COV_N_est(k, k2, n);
                end
            end
        end
    end
    Sigma_N_est = COV_N_est / p_Kotz;
    
    %% Run algorithm
    cost_function = [];
    cost_const = N * (K / 2 * log(pi) + gammaln(K/2/beta) + K / 2 / beta * log(2) - log(beta) - gammaln(K/2));
    for iter = 1:Params.maxIter
    
        cost = cost_const;
        for n = 1:N
            Yn = squeeze(Y(n, :, :)); % TxK
            Zn = Yn / Sigma_N_est(:, :, n); % TxK
            cost = cost + .5 * logdet(Sigma_N_est(:, :, n)) + 0.5 * mean(dot(Zn, Yn, 2).^beta);
        end
        for k = 1:K
            cost = cost - log(abs(det(W(:, :, k))));
        end
        cost_function = [cost_function, cost];
        if length(cost_function) > 1 && cost_function(end) > cost_function(end-1)
            stepsize = max(stepsize*Params.stepsizeDecrease, Params.minStepsize);
        end
    
        W_old = W;
        Q = 0;
        R = 0;
        for n = 1:N
            [Dn, Q, R] = decouple_trick(W, n, Q, R); % N*K
    
            Yn = squeeze(Y(n, :, :)); % TxK
            Zn = Yn / Sigma_N_est(:, :, n); % TxK
            ySiy = dot(Zn, Yn, 2).^(beta - 1); % Tx1
            for k = 1:K
                phi_nk = Zn(:, k) .* ySiy; % Tx1
                dWnk = (beta/T) * X(:, :, k) * phi_nk - Dn(:, k) / (W(n, :, k) * Dn(:, k)); % Nx1
    
                dWnk = dWnk - W(n, :, k) * dWnk * W(n, :, k)'; % Nx1
                dWnk = dWnk / norm(dWnk); % check Li's 2010 TSP paper
                temp = W(n, :, k) - stepsize * dWnk';
                W(n, :, k) = temp / norm(temp); % 1xN
    
                % update variables that depend on Wnk
                Yn(:, k) = W(n, :, k) * X(:, :, k);
                for k2 = 1:K
                    Sigma_N_est(k, k2, n) = W(n, :, k) * Rx{k, k2} * W(n, :, k2)' / p_Kotz;
                    Sigma_N_est(k2, k, n) = Sigma_N_est(k, k2, n);
                end
                if Params.updateCOV && k < K
                    Zn = Yn / Sigma_N_est(:, :, n); % TxK
                    ySiy = dot(Zn, Yn, 2).^(beta - 1); % Tx1
                end
            end % k
    
            Y(n, :, :) = Yn;
        end % n
    
        currentChange = 0;
        for k = 1:K
            currentChange = max(currentChange, max(1-abs(diag(W_old(:, :, k)*W(:, :, k)'))));
        end
    
        if currentChange < Params.termThreshold && iter > Params.minIter
            break;
        end
    
        if Params.checkpoints && mod(iter, 25) == 0
            save(['cIVA_', num2str(iter), '_', char(datetime('now'), 'MM_dd_yyyy_HH_mm_ss'), '.mat'], 'W', 'cost_function', ...
                'stepsize', 'currentChange', 'COV_N_est', 'dWnk')
        end
    end
    
    %% Dewhitening
    if Params.whiten
        for k = 1:K
            W(:, :, k) = W(:, :, k) * V(:, :, k);
        end
    end
    
    % Resort order of SCVs from most to least ill-conditioned
    if Params.sortSCVs
        det_COV = zeros(N, 1); 
        for n = 1:N
            det_COV(n) = det(COV_N_est(:, :, n));
        end
        [~, isort] = sort(det_COV);
        % COV_N_est=COV_N_est(:,:,isort);
        for k = 1:K
            W(:, :, k) = W(isort, :, k);
        end
    end

%% Return
return;
% END OF iva

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

function v = logdet(A, op)
%LOGDET Computation of logarithm of determinant of a matrix
%
%   v = logdet(A);
%       computes the logarithm of determinant of A. 
%
%       Here, A should be a square matrix of double or single class.
%       If A is singular, it will returns -inf.
%
%       Theoretically, this function should be functionally 
%       equivalent to log(det(A)). However, it avoids the 
%       overflow/underflow problems that are likely to 
%       happen when applying det to large matrices.
%
%       The key idea is based on the mathematical fact that
%       the determinant of a triangular matrix equals the
%       product of its diagonal elements. Hence, the matrix's
%       log-determinant is equal to the sum of their logarithm
%       values. By keeping all computations in log-scale, the
%       problem of underflow/overflow caused by product of 
%       many numbers can be effectively circumvented.
%
%       The implementation is based on LU factorization.
%
%   v = logdet(A, 'chol');
%       If A is positive definite, you can tell the function 
%       to use Cholesky factorization to accomplish the task 
%       using this syntax, which is substantially more efficient
%       for positive definite matrix. 
%
%   Remarks
%   -------
%       logarithm of determinant of a matrix widely occurs in the 
%       context of multivariate statistics. The log-pdf, entropy, 
%       and divergence of Gaussian distribution typically comprises 
%       a term in form of log-determinant. This function might be 
%       useful there, especially in a high-dimensional space.       
%
%       Theoretially, LU, QR can both do the job. However, LU 
%       factorization is substantially faster. So, for generic
%       matrix, LU factorization is adopted. 
%
%       For positive definite matrices, such as covariance matrices,
%       Cholesky factorization is typically more efficient. And it
%       is STRONGLY RECOMMENDED that you use the chol (2nd syntax above) 
%       when you are sure that you are dealing with a positive definite
%       matrix.
%
%   Examples
%   --------
%       % compute the log-determinant of a generic matrix
%       A = rand(1000);
%       v = logdet(A);
%
%       % compute the log-determinant of a positive-definite matrix
%       A = rand(1000);
%       C = A * A';     % this makes C positive definite
%       v = logdet(C, 'chol');
%

%   Copyright 2008, Dahua Lin, MIT
%   Email: dhlin@mit.edu
%
%   This file can be freely modified or distributed for any kind of 
%   purposes.
%

%% argument checking

assert(isfloat(A) && ndims(A) == 2 && size(A,1) == size(A,2), ...
    'logdet:invalidarg', ...
    'A should be a square matrix of double or single class.');

if nargin < 2
    use_chol = 0;
else
    assert(strcmpi(op, 'chol'), ...
        'logdet:invalidarg', ...
        'The second argument can only be a string ''chol'' if it is specified.');
    use_chol = 1;
end

%% computation

if use_chol
    v = 2 * sum(log(diag(chol(A))));
else
    [L, U, P] = lu(A);
    du = diag(U);
    c = det(P) * prod(sign(du));
    v = log(c) + sum(log(abs(du)));
end
