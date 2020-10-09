function [W,cost,Sigma_N,isi] = icatb_iva_second_order(X,varargin)
% [W,cost,Sigma_N,isi] = iva_second_order(X,varargin)
%
% Implementations of all the second-order (Gaussian) independent vector
% analysis (IVA) algorithms.  Namely real-valued and complex-valued with
% circular and non-circular using gradient, Newton, and quasi-Newton
% optimizations.
%
% Input:
% X - data observations from K data sets, i.e. X{k}=A{k}*S{k}, where A{k}
% is an N x N unknowned invertible mixing matrix and S{k} is N x T  matrix
% with the nth row corresponding to T samples of the nth source in the kth
% dataset.  For IVA it is assumed that the source is statistically
% independent of all the sources within a dataset and exactly dependent on
% at most one source in each of the other datasets. The data, X, can be
% either a cell array as indicated or a 3-dimensional matrix of dimensions
% N x T x K.  The latter enforces the assumption of an equal number of
% samples in each dataset.
%
% Output:
% W - the estimated demixing matrices so that ideally W{k}*A{k} = P*D{k}
% where P is any arbitrary permutation matrix and D{k} is any diagonal
% invertible (scaling) matrix.  Note P is common to all datasets; this is
% to indicate that the local permuation ambiguity between dependent sources
% across datasets should ideally be resolved by IVA.
%
% cost - the cost for each iteration
%
% isi - joint inter-symbol-interference is available if user supplies true
% mixing matrices for computing a performance metric
%
% Optional input pairs and default values are given in list below:
%
%    'opt_approach','newton', ... % optimization type: gradient, (newton), quasi
%    'complex_valued',false, ... % if any input is complex or quasi approach used then setting is forced to true
%    'circular',false, ... % set to true to only consider circular for complex-valued cost function
%    'whiten',true, ... % whitening is optional (except for quasi approach it is required)
%    'verbose',false, ... % verbose true enables print statements
%    'A',[], ... % true mixing matrices A, automatically sets verbose
%    'W_init',[], ... % initial estimates for demixing matrices in W
%    'jdiag_initW',false, ... % use CCA (K=2) or joint diagonalization (K>2)
%    'maxIter',512, ... % max number of iterations
%    'WDiffStop',1e-6, ... % stopping criterion
%    'alpha0',1.0 ... % initial step size scaling (will be doubled for complex-valued)
%
% Example call:
% W=iva_second_order(X,'opt_approach','newton','circular',true)
%
% Given that X is complex-valued this uses Newton optimization for circular
% Gaussian IVA cost function.  If the last argument is set to false then
% the cost-function is the more general non-circular Gaussian IVA cost
% function.
%
% To test code call function with no input arguments: iva_second_order
%
% Coded by Matthew Anderson (matt.anderson at umbc.edu)
%
% References:
%
% [1] M. Anderson, X.-L. Li, & T. Adal??, "Nonorthogonal Independent Vector Analysis Using Multivariate Gaussian Model," LNCS: Independent Component Analysis and Blind Signal Separation, Latent Variable Analysis and Signal Separation, Springer Berlin / Heidelberg, 2010, 6365, 354-361
% [2] M. Anderson, T. Adal??, & X.-L. Li, "Joint Blind Source Separation of Multivariate Gaussian Sources: Algorithms and Performance Analysis," IEEE Trans. Signal Process., 2012, 60, 1672-1683
% [3] M. Anderson, X.-L. Li, & T. Adal??, "Complex-valued Independent Vector Analysis: Application to Multivariate Gaussian Model," Signal Process., 2012, 1821-1831

% Version 01 - 20120913 - Initial publication
% Version 02 - 20120919 - Using decouple_trick function, bss_isi, whiten, &
%                         cca subfunctions

global JBSS_SOS_EXIST

if isempty(JBSS_SOS_EXIST)
    JBSS_SOS_EXIST=exist('jbss_sos.m','file');
end
if nargin==0
    help iva_second_order
    test_iva_second_order
    return
end


%% Gather Options

% build default options structure
options=struct( ...
    'opt_approach','newton', ... % optimization type: gradient, (newton), quasi
    'complex_valued',false, ... % if any input is complex or quasi approach used then setting is forced to true
    'circular',false, ... % set to true to only consider circular for complex-valued cost function
    'whiten',true, ... % whitening is optional (except for quasi approach it is required)
    'verbose',false, ... % verbose true enables print statements
    'A',[], ... % true mixing matrices A, automatically sets verbose
    'W_init',[], ... % initial estimates for demixing matrices in W
    'jdiag_initW',false, ... % use CCA (K=2) or joint diagonalization (K>2)
    'maxIter',512, ... % max number of iterations
    'WDiffStop',1e-6, ... % stopping criterion
    'alpha0',1.0 ... % initial step size scaling (will be doubled for complex-valued)
    );

% load in user supplied options
options=getopt(options,varargin{:});

supplyA=~isempty(options.A); % set to true if user has supplied true mixing matrices
blowup = 1e3;
alphaScale=0.9; % alpha0 to alpha0*alphaScale when cost does not decrease
alphaMin=options.WDiffStop; % alpha0 will max(alphaMin,alphaScale*alpha0)
outputISI=false;
opt_approach=find(strcmp(options.opt_approach, ...
    {'gradient','newton','quasi'}),1); % 1='gradient', 2='newton', 3='quasi'
options.whiten=(options.whiten || (opt_approach==3)); % whitening required for real-valued quasi & complex-valued gradient approach


%% Create cell versions
% Faster than using the more natural 3-dimensional version
if iscell(options.A)
    K=length(options.A);
    A=zeros([size(options.A{1}) K]);
    for k=1:K
        A(:,:,k)=options.A{k};
    end
    options.A=A;
    options.complex_valued = options.complex_valued || any(imag(A(:)));
    clear A
end

if ~iscell(X)
    options.complex_valued = options.complex_valued || any(imag(X(:)));
    options.whiten = options.whiten || (options.complex_valued &&  opt_approach==1);
    if options.whiten
        [X,V]=whiten(X);
    end
    [N,T,K]=size(X); % the input format insures all datasets have equal number of samples
    %     Rx=cell(K,K);
    %     for k1=1:K
    %         Rx{k1,k1}=1/T*(X(:,:,k1)*X(:,:,k1)');
    %         for k2=(k1+1):K
    %             Rx{k1,k2}=1/T*(X(:,:,k1)*X(:,:,k2)');
    %             Rx{k2,k1}=Rx{k1,k2}';
    %         end
    %     end
    
    Rx = computeCrossCov(X); % Optimize cross covariance calculation
    
    if options.complex_valued && ~options.circular
        % complex-valued non-circular cost
        % compute data pseudo cross-covariance matrices
        Px=cell(K,K);
        for k1=1:K
            Px{k1,k1}=1/T*(X(:,:,k1)*X(:,:,k1).');
            for k2=(k1+1):K
                Px{k1,k2}=1/T*(X(:,:,k1)*X(:,:,k2).');
                Px{k2,k1}=Px{k1,k2}.';
            end % k2
        end % k1
    end
else % iscell
    [rowX,colX]=size(X);
    
    % then this is cell array with each entry being X{k} being a dataset
    % for now assume N is fixed
    X=X(:);
    
    K=max(colX,rowX);
    [N,T]=size(X{1});
    if ~options.complex_valued
        for k=1:K
            if T~=size(X{k},2)
                error('Each dataset must have the same number of samples.')
            end
            options.complex_valued = options.complex_valued || any(imag(X{k}(:)));
            if options.complex_valued
                break
            end
        end % k
    end
    
    options.whiten = options.whiten || (options.complex_valued &&  opt_approach==1);
    
    if options.whiten
        [X,V_cell]=whiten(X);
        V=zeros(N,N,K);
        for k=1:K
            V(:,:,k)=V_cell{k};
        end
        clear V_cell
    end
    Rx=cell(K,K);
    for k1=1:K
        
        options.complex_valued = options.complex_valued || any(imag(X{k1}(:)));
        Rx{k1,k1}=1/T*(X{k1}*X{k1}');
        for k2=(k1+1):K
            Rx{k1,k2}=1/T*(X{k1}*X{k2}');
            Rx{k2,k1}=Rx{k1,k2}';
        end
    end
    if options.complex_valued && ~options.circular
        % complex-valued non-circular cost
        % compute data pseudo cross-covariance matrices
        Px=cell(K,K);
        for k1=1:K
            Px{k1,k1}=1/T*(X{k1}*X{k1}.');
            for k2=(k1+1):K
                Px{k1,k2}=1/T*(X{k1}*X{k2}.');
                Px{k2,k1}=Px{k1,k2}.';
            end % k2
        end % k1
    end
end

%% Initialize W
if ~isempty(options.W_init)
    W=options.W_init;
    if ~options.complex_valued && any(imag(W(:))) && opt_approach==1
        warning(['Supplied initial W is complex valued and ' ...
            'gradient option approach is set; ' ...
            'need to set whitening option to true for best performance.']) %#ok<WNTAG>
    end
    if size(W,3)==1 && size(W,1)==N && size(W,2)==N
        W=repmat(W,[1,1,K]);
    end
    if options.whiten
        for k=1:K
            W(:,:,k)=W(:,:,k)/V(:,:,k);
        end
    end
else
    if options.jdiag_initW
        
        if K>2
            % initialize with multi-set diagonalization (orthogonal solution)
            % Assume N is same for all K datasets
            if ~JBSS_SOS_EXIST
                error('To use this option requires another package.')
            end
            W=jbss_sos(X,0,'whole');
        else % if K==2 then CCA is sufficient (should be few to no iterations required)
            W_cell=cca(X);
            W=zeros(N,N,K);
            for k=1:K
                W(:,:,k)=W_cell{k};
            end
        end
    else % randomly initialize
        W=randn(N,N,K);
        if options.complex_valued
            W=W+1j*randn(N,N,K);
        end
        for k = 1:K
            W(:,:,k)=(sqrtm(W(:,:,k)*W(:,:,k)'))\W(:,:,k);
        end
    end
end

%% When mixing matrix A supplied
% verbose is set to true
% outputISI can be computed if requested
% A matrix is conditioned by V if data is whitened
if supplyA
    % only reason to supply A matrices is to display running performance
    options.verbose=true;
    if nargout>3
        outputISI=true;
        isi=zeros(1,options.maxIter);
    end
    if options.whiten
        for k = 1:K
            options.A(:,:,k) = V(:,:,k)*options.A(:,:,k);
        end
    end
end

%% Check rank of data-covariance matrix
% should be full rank, if not we inflate
% this is ad hoc and no robustness is assured by diagonal loading procedure
% Rxall=cell2mat(Rx);
% r=rank(Rxall);
% if r<(N*K)
%     % inflate Rx
%     eigRx=svd(Rxall);
%     Rxall=Rxall+eigRx(r)*eye(N*K); %diag(linspace(0.1,0.9,N*K)); % diagonal loading
%     Rx=mat2cell(Rxall,N*ones(K,1),N*ones(K,1));
% end
%
% clear Rxall;

%% Initialize some local variables
cost=zeros(1,options.maxIter);
cost_const=K*log(2*pi*exp(1)); % local constant

%% Initializations based on real vs complex-valued
if options.complex_valued
    grad=complex(zeros(N,K));
    if opt_approach==2
        HA=zeros(N*K,N*K);
    end
    
    % Double step-size for complex-valued gradient optimization
    if opt_approach==1
        options.alpha0 = 2*options.alpha0;
    end
else % real-valued
    grad=zeros(N,K);
    if opt_approach==2
        H=zeros(N*K,N*K);
    end
end


isCellX = iscell(X);

%clear X;



%% Main Iteration Loop
for iter = 1:options.maxIter
    termCriterion=0;
    
    % Some additional computations of performance via ISI when true A is supplied
    if supplyA
        [amari_avg_isi,amari_joint_isi]=bss_isi(W,options.A);
        if outputISI
            isi(iter)=amari_joint_isi;
        end
    end
    
    W_old=W; % save current W as W_old
    cost(iter)=0;
    for k=1:K
        cost(iter)=cost(iter)-log(abs(det(W(:,:,k))));
    end
    
    Q=0; R=0;
    %% Loop over each SCV
    for n=1:N
        Wn=W(n,:,:);
        Wn=conj(Wn(:));
        
        % Efficient version of Sigma_n=Yn*Yn'/T;
        Sigma_n=eye(K); %
        for k1=1:K
            for k2=k1:K
                Sigma_n(k1,k2)=W(n,:,k1)*Rx{k1,k2}*W(n,:,k2)';
                if k1~=k2
                    Sigma_n(k2,k1)=Sigma_n(k1,k2)';
                end
            end % k2
        end %k1
                
        
        if options.complex_valued && ~options.circular
            % complex-valued non-circular cost
            % compute SCV pseudo cross-covariance matrices
            SigmaP_n=zeros(K); % pseudo = Yn*Yn.'/T;
            for k1=1:K
                for k2=k1:K
                    SigmaP_n(k1,k2)=W(n,:,k1)*Px{k1,k2}*W(n,:,k2).';
                    if k1~=k2
                        SigmaP_n(k2,k1)=SigmaP_n(k1,k2);
                    end
                end % k2
            end %k1
            Sigma_Total=[Sigma_n SigmaP_n; SigmaP_n' conj(Sigma_n)];
            cost(iter)=cost(iter)+0.5*log(det(Sigma_Total));
        else
            cost(iter)=cost(iter)+0.5*(cost_const+log(det(Sigma_n)));
        end
        
        if options.complex_valued && ~options.circular
            inv_Sigma_n=inv(Sigma_Total);
            P=inv_Sigma_n(1:K,1:K);
            Pt=inv_Sigma_n(1:K,K+(1:K));
        else
            inv_Sigma_n=inv(Sigma_n);
        end
        [hnk,Q,R]=decouple_trick(W,n,Q,R);
        
        for k=1:K
            % Analytic derivative of cost function with respect to vn
            % Code below is efficient implementation of computing the gradient, that is
            % independent of T
            grad(:,k)=-hnk(:,k)/(W(n,:,k)*hnk(:,k));
            if options.complex_valued && ~options.circular
                for kk=1:K
                    grad(:,k)=grad(:,k)+ ...
                        Rx{k,kk}*W(n,:,kk)'*P(k,kk)'+Px{k,kk}*W(n,:,kk).'*Pt(k,kk)';
                end
            else
                
                for kk=1:K
                    grad(:,k)=grad(:,k)+ ...
                        Rx{k,kk}*W(n,:,kk)'*inv_Sigma_n(kk,k);
                end
                
            end
            
            if opt_approach==1 || ... % real/complex gradient
                    (opt_approach==3 && ~options.complex_valued) % real quasi newton
                wnk=W(n,:,k)';
                if opt_approach==1
                    gradNorm=vecnorm(grad(:,k)); % normalize gradient
                    gradNormProj=vecnorm(gradNorm-wnk'*gradNorm*wnk); % non-colinear direction normalized
                    W(n,:,k)=vecnorm(wnk-options.alpha0*gradNormProj)';
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%% The following computations are potentially unnecessary
                    %%%%% see that the quasi-newton approach does not do this
                    %%%%% step.
                    for kk=1:K
                        Sigma_n(k,kk)=W(n,:,k)*Rx{k,kk}*W(n,:,kk)'; % = Yn*Yn'/T;
                    end
                    Sigma_n(:,k)=Sigma_n(k,:)';
                    if options.complex_valued && ~options.circular
                        for kk=1:K
                            SigmaP_n(k,kk)=W(n,:,k)*Px{k,kk}*W(n,:,kk).'; % = Yn*Yn.'/T
                        end
                        SigmaP_n(:,k)=SigmaP_n(k,:).';
                        Sigma_Total=[Sigma_n SigmaP_n; SigmaP_n' conj(Sigma_n)];
                        inv_Sigma_n=inv(Sigma_Total);
                    else
                        inv_Sigma_n=inv(Sigma_n);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                elseif opt_approach==3 % real-valued quasi-newton (block newton)
                    % Hessian inverse computed using matrix inversion lemma
                    Hinv=1/inv_Sigma_n(k,k)*(eye(N) - ...
                        hnk(:,k)*hnk(:,k)'/(inv_Sigma_n(k,k)+1/(hnk(:,k)'*wnk)^2));
                    
                    % Hinv = bsxfun(@minus, eye(N), hnk(:,k)*hnk(:,k)'/(inv_Sigma_n(k,k)+1/(hnk(:,k)'*wnk)^2))/inv_Sigma_n(k,k);
                    
                    % Block-Newton update of Wnk
                    W(n,:,k)=vecnorm(wnk-options.alpha0*(Hinv*grad(:,k)))';
                end
            end
        end % k
        
        if opt_approach==2 || (opt_approach==3 && options.complex_valued)
            %% Compute SCV Hessian
            for k1=1:K
                if options.complex_valued
                    if opt_approach==2 % complex-valued Newton
                        HA((k1-1)*N+(1:N),(k1-1)*N+(1:N)) = ...
                            conj(hnk(:,k1))*hnk(:,k1)'/(hnk(:,k1)'*W(n,:,k1)')^2;
                        if ~options.circular % complex-valued Newton non-circular
                            HA((k1-1)*N+(1:N),(k1-1)*N+(1:N)) = ...
                                HA((k1-1)*N+(1:N),(k1-1)*N+(1:N)) + ...
                                Pt(k1,k1)*conj(Px{k1,k1});
                        end
                    end
                    H((k1-1)*N+(1:N),(k1-1)*N+(1:N)) = ...
                        conj(inv_Sigma_n(k1,k1)*Rx{k1,k1});
                else % real-valued Newton
                    H((k1-1)*N+(1:N),(k1-1)*N+(1:N)) = ...
                        inv_Sigma_n(k1,k1)*Rx{k1,k1} + ...
                        hnk(:,k1)*hnk(:,k1)'/(hnk(:,k1)'*W(n,:,k1)')^2;
                end
                for k2=(k1+1):K
                    if opt_approach==2 && options.complex_valued && ~options.circular
                        % complex-valued Newton non-circular
                        HA((k1-1)*N+(1:N),(k2-1)*N+(1:N)) = ...
                            Pt(k2,k1)*conj(Px{k1,k2});
                        HA((k2-1)*N+(1:N),(k1-1)*N+(1:N)) = ...
                            Pt(k1,k2)*conj(Px{k2,k1});
                        Hs=conj(P(k2,k1)*Rx{k1,k2});
                    else
                        Hs=inv_Sigma_n(k1,k2)*Rx{k2,k1}.';
                    end
                    H((k1-1)*N+(1:N),(k2-1)*N+(1:N))=Hs;
                    H((k2-1)*N+(1:N),(k1-1)*N+(1:N))=Hs';
                end % k2
            end % k1
            
            %% Alternative approach to calculating H
            % H=kron(inv_Sigma_n,ones(N)).*Rxall;
            % for k=1:K
            %    H((k-1)*N+(1:N),(k-1)*N+(1:N))=H((k-1)*N+(1:N),(k-1)*N+(1:N)) + ...
            %         hnk(:,k)*hnk(:,k)'/(hnk(:,k)'*W(n,:,k)')^2;
            % end
            
            %% Newton Update
            if ~options.complex_valued % real-valued Newton               
                ysol = cholsol(H, grad(:));            
                %[ysol, flag] = gmres(H, grad(:), [], 1e-6, 50);
                %if (flag ~= 0)
                %    ysol = (H\grad(:));
                %end
                
                Wn=Wn-options.alpha0*ysol; % much faster inv(H)*grad(:)
                %[ysol, flag] = gmres(H, grad(:), [], 1e-6, 20);
                %Wn=Wn-options.alpha0*ysol; % much faster inv(H)*grad(:)
            else % complex-valued
                if opt_approach == 2 % complex-valued Newton
                    % improve stability by using 0.99 in Hessian
                    Wn=Wn-options.alpha0*(conj(H)-0.99*conj(HA)*(H\HA))\(grad(:)-conj(HA)*(H\conj(grad(:))));
                elseif opt_approach == 3 % complex-valued Quasi-Newton
                    Wn=Wn-options.alpha0*conj(H)\grad(:);
                else
                    error('Should not happen.')
                end
                Wn=conj(Wn);
            end
            
            %% Store Updated W
            Wn=(reshape(Wn,N,K));
            for k=1:K
                W(n,:,k)=vecnorm(Wn(:,k)).';
            end % k
        end % opt_approach==2 || (opt_approach==3 && options.complex_valued)
    end %n
    
    %    for k=1:K
    %       termCriterion = max(termCriterion,max(1-abs(diag(W_old(:,:,k)*W(:,:,k)'))));
    %    end % k
    
    % this is what we would like to change it to.
    for k=1:K
        termCriterion = max(termCriterion,norm(W(:,:,k) - W_old(:,:,k),'fro'));
    end % k
    
    
    %% Decrease step size alpha if cost increased from last iteration
    if iter>1
        if cost(iter)>cost(iter-1)
            options.alpha0=max(alphaMin,alphaScale*options.alpha0);
        end
    end
    
    
    if (iter==1)
        costCriterion = 1;
    else
        costCriterion = abs(cost(iter-1)-cost(iter))/abs(cost(iter));
    end
    
    %% Check the termination condition
    if (termCriterion < options.WDiffStop || iter == options.maxIter || (costCriterion < options.WDiffStop))
        break;
    elseif termCriterion > blowup || isnan(cost(iter))
        for k = 1:K
            W(:,:,k) = eye(N) + 0.1*randn(N);
        end
        if options.verbose
            fprintf('\n W blowup, restart with new initial value.');
        end
    end
    
    %% Display Iteration Information
    if options.verbose
        if supplyA
            fprintf('\n Step %d: W change: %f, Cost: %f, Avg ISI: %f, Joint ISI: %f',  ...
                iter, termCriterion,cost(iter),amari_avg_isi,amari_joint_isi);
        else
            fprintf('\n Step %d: W change: %f, Cost: %f', iter, termCriterion,cost(iter));
        end
    end % options.verbose
end % iter

%% Finish Display
if iter==1 && options.verbose
    if supplyA
        fprintf('\n Step %d: W change: %f, Cost: %f, Avg ISI: %f, Joint ISI: %f',  ...
            iter, termCriterion,cost(iter),amari_avg_isi,amari_joint_isi);
    else
        fprintf('\n Step %d: W change: %f, Cost: %f', iter, termCriterion,cost(iter));
    end
end % options.verbose

if options.verbose
    fprintf('\n');
end

%% Clean-up Outputs
cost=cost(1:iter);
if outputISI
    isi=isi(1:iter);
end

if options.whiten
    for k=1:K
        W(:,:,k)=W(:,:,k)*V(:,:,k);
    end
else % no prewhitening
    %% Scale demixing vectors to generate unit variance sources
    for n=1:N
        for k=1:K
            W(n,:,k)=W(n,:,k)/sqrt(W(n,:,k)*Rx{k,k}*W(n,:,k)');
        end
    end
end

%% Resort order of SCVs
% Order the components from most to least ill-conditioned

% First, compute the data covariance matrices (by undoing any whitening)
if options.whiten
    for k1=1:K
        Rx{k1,k1}=(V(:,:,k1)\Rx{k1,k1})/V(:,:,k1)';
        for k2=(k1+1):K
            Rx{k1,k2}=(V(:,:,k1)\Rx{k1,k2})/V(:,:,k2)';
            Rx{k2,k1}=Rx{k1,k2}';
        end % k2
    end % k1
end

% Second, compute the determinant of the SCVs
detSCV=zeros(1,N);
Sigma_N=zeros(K,K,N);
for n=1:N
    % Efficient version of Sigma_n=Yn*Yn'/T;
    Sigma_n=zeros(K); %
    for k1=1:K
        for k2=k1:K
            Sigma_n(k1,k2)=W(n,:,k1)*Rx{k1,k2}*W(n,:,k2)';
            Sigma_n(k2,k1)=Sigma_n(k1,k2)';
        end % k2
    end %k3
    Sigma_N(:,:,n)=Sigma_n;
    
    if any(imag(Sigma_n(:))~=0) && ~options.circular
        SigmaP_n=zeros(K); % pseudo = Yn*Yn.'/T;
        for k1=1:K
            for k2=k1:K
                SigmaP_n(k1,k2)=W(n,:,k1)*Px{k1,k2}*W(n,:,k2).';
                if k1~=k2
                    SigmaP_n(k2,k1)=SigmaP_n(k1,k2);
                end
            end % k2
        end %k1
        detSCV(n)=det([Sigma_n SigmaP_n; SigmaP_n' conj(Sigma_n)]);
    else
        detSCV(n)=det(Sigma_n);
    end
end

% Third, sort and apply
[ddddd,isort]=sort(detSCV);
Sigma_N=Sigma_N(:,:,isort);
for k=1:K
    W(:,:,k)=W(isort,:,k);
end

%% Cast output
if isCellX %iscell(X)
    Wmat=W;
    W=cell(K,1);
    for k=1:K
        W{k}=Wmat(:,:,k);
    end
end

%% Return
return;
% END OF IVA_MG_OPTS

function test_iva_second_order
%% built-in test function for iva_second_order
N = 10;  % number of sources
K = 10;  % number of groups
T = 1000;   % sample size

% generate the mixtures
S = zeros(N,T,K);
for n = 1 : N
    
    temp1 = randn(K,T);
    temp  = zeros(K,T);
    B = randn(K,K,3);
    for p = 0 : 2
        for t = 3 : T
            temp(:,t) = temp(:,t) + B(:,:,p+1)*temp1(:,t-p); % introduce nonwhiteness and spatial correlation
        end
    end
    
    for k = 1 : K
        S(n,:,k) = temp(k,:);
        S(n,:,k) = S(n,:,k) - mean(S(n,:,k));
        S(n,:,k) = S(n,:,k)/sqrt(var(S(n,:,k)));
    end
end
A = randn(N,N,K);
x = zeros(N,T,K);
for k = 1 : K
    x(:,:,k) = A(:,:,k)*S(:,:,k);
end


% separation
figure;

W=iva_second_order(x,'jdiag_initW',false,'opt_approach','newton');
% show results
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
[ddddd,imax]=max(T1);
ind=sub2ind([N N],1:N,imax);
P(ind)=1;
T1=P*T1;
imagesc(1:N,1:N,T1)
caxis([0 1])
axis xy
colormap('bone')
title('joint global matrix')
colorbar
disp('Ideally image is identity matrix.')


%% Available options:
%   'opt_approach','newton', ... % optimization type: gradient, (newton), quasi
%    'verbose',false, ... % verbose true enables print statements
%    'A',[], ... % true mixing matrices A, automatically sets verbose
%    'W_init',[], ... % initial estimates for demixing matrices in W
%    'jdiag_initW',true, ... % use CCA (K=2) or joint diagonalization (K>2)
%    'maxIter',512, ... % max number of iterations
%    'WDiffStop',1e-6, ... % stopping criterion
%    'alpha0',1.0 ... % initial step size scaling

% %% Example calls to iva_second_order with various options
% % note that the inputs could be made into a cell array with K entries
% W=iva_second_order(X,'jdiag_initW',false,'opt_approach','newton');
% [~,isi]=bss_isi(W,A,S);
% disp(['Newton IVA-G ISI: ' num2str(isi)])
W=iva_second_order(x,'A',A,'jdiag_initW',false,'opt_approach','gradient','whiten',true);
[ddddd,isi]=bss_isi(W,A,S);
% disp(['Block-Newton IVA-G ISI: ' num2str(isi)])
% X_cell=cell(K,1);
% for k=1:K
%    X_cell{k}=X(:,:,k);
% end
% W_cell=iva_second_order(X_cell,'jdiag_initW',false,'opt_approach','newton');
% for k=1:K
%    W(:,:,k)=W_cell{k};
% end
% [~,isi]=bss_isi(W,A,S);
% disp(['Block-Newton IVA-G ISI: ' num2str(isi)])
%

% %% Generate test samples:
% S=complex(zeros(N,T,K)); % X(:,:,k)=A(:,:,k)*S(:,:,k);
%
% for n=1:N
%    Z = zeros(K,T);
%    B=randn(K);
%    B=B*B';
%
%    [V,D]=eig(B);
%    C=V*sqrt(D);
%    absrho=sqrt(rand(1,K));
%    angrho=rand(1,K)*2*pi;
%    rho=absrho.*exp(1j*angrho);
%    for k=1:K
%       Cggd=[1 rho(k); conj(rho(k)) 1];
%       [X]=randggd_complex(T,Cggd,1);
%       Z(k,:)=X(1,:);
%    end
%    Z=C*Z;
%    %    if n==1
%    %       for k=2:K
%    %          Z(k,:)=Z(1,:);
%    %       end
%    %    end
%    S(n,:,:)=Z.';
% end
% A=randn(N,N,K)+1j*randn(N,N,K);
% X=S;
% for k=1:K
%    X(:,:,k)=A(:,:,k)*S(:,:,k);
% end
%
%
% %% Available options:
% %   'opt_approach','newton', ... % optimization type: gradient, (newton), quasi
% %    'verbose',false, ... % verbose true enables print statements
% %    'A',[], ... % true mixing matrices A, automatically sets verbose
% %    'W_init',[], ... % initial estimates for demixing matrices in W
% %    'jdiag_initW',true, ... % use CCA (K=2) or joint diagonalization (K>2)
% %    'maxIter',512, ... % max number of iterations
% %    'WDiffStop',1e-6, ... % stopping criterion
% %    'alpha0',1.0 ... % initial step size scaling
%
% %% Example calls to iva_second_order with various options
% % note that the inputs could be made into a cell array with K entries
% W=iva_second_order(X,'A',A,'jdiag_initW',false,'opt_approach','newton','circular',true);
% [~,isi]=bss_isi(W,A,S);
% disp(['Circular IVA-G ISI: ' num2str(isi)])
% W=iva_second_order(X,'A',A,'jdiag_initW',false,'opt_approach','newton');
% [~,isi]=bss_isi(W,A,S);
% disp(['Non-Circ IVA-G ISI: ' num2str(isi)])

%%
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
    % mag=sqrt(sum(vec.*conj(vec)));
    mag = sqrt(sum(bsxfun(@times, vec, conj(vec))));
else
    mag=vec.*conj(vec);
    for ii=1:length(varargin)
        mag=mag+varargin{ii}.*conj(varargin{ii});
    end
    mag=sqrt(mag);
end
%return

function [uvec,mag]=vecnorm(vec)
% [vec,mag]=vecnorm(vec)
% Returns the vector normalized by 2-norm or magnitude of vector.
% vec has size n by m represents m vectors of length n (i.e. m
% column-vectors).
[n,m]=size(vec);
if n==1
    disp('vecnorm operates on column vectors, input appears to have dimension of 1')
end

%uvec=zeros(n,m);
mag=vecmag(vec); % returns a 1 x m row vector
% for ii=1:size(vec,1)
%     uvec(ii,:)=vec(ii,:)./mag;
% end

uvec = bsxfun(@rdivide, vec, mag);

% Equivalent to: uvec=vec./repmat(mag,size(vec,1),1);

% Which implementation is optimal depends on optimality criterion (memory
% vs. speed), this version uses the former criterion.
%return


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
    
    %         if (n == 1)
    %
    %             parfor k = 1:K
    %                 %iQ = invQ;
    %                 Wtilde = W;
    %                 Wtilde=Wtilde(2:N,:,k);
    %                 [tmpQ,tmpR]=qr(Wtilde');
    %                 invQ(:, :, k) = tmpQ;
    %                 R(:,:,k) = tmpR;
    %                 h(:,k)=tmpQ(:,end);
    %             end
    %
    %         else
    %
    %             parfor k = 1:K
    %                 WT = W;
    %                 n_last=n-1;
    %                 e_last = zeros(N-1,1);
    %                 e_last(n_last) = 1;
    %                 [invQ(:,:,k),R(:,:,k)]=qrupdate(invQ(:,:,k),R(:,:,k),-WT(n,:,k)',e_last);
    %                 [invQ(:,:,k),R(:,:,k)]=qrupdate(invQ(:,:,k),R(:,:,k),WT(n_last,:,k)',e_last);
    %                 tmp = invQ(:, :, k);
    %                 h(:,k)=tmp(:,end);
    %             end
    %
    %         end
    
    
    for k=1:K
        if n==1
            Wtilde=W(2:N,:,k);
            [invQ(:,:,k),R(:,:,k)]=qr(Wtilde');
        else
            n_last=n-1;
            e_last = zeros(N-1,1);
            e_last(n_last) = 1;
            [invQ(:,:,k),R(:,:,k)]=qrupdate(invQ(:,:,k),R(:,:,k), (W(n_last,:,k)-W(n,:,k))',e_last);
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
        [Q,dddd]=qr(W([(1:n-1) (n+1:N)],:,k)');
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
                [dddd,colMaxG]=max(Gabs,[],2);
                if length(unique(colMaxG))~=Np
                    % solution is failure in strictest sense
                    success=false;
                end
            else
                [dddd,colMaxG_k]=max(Gabs,[],2);
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
            [dddd,colMaxG]=max(Gabs(1:Nk,1:Nk),[],2);
            if length(unique(colMaxG))~=Nk
                % solution is a failure in a strict sense
                success=false;
            end
        elseif success
            if nargin>=4
                [dddd,colMaxG_k]=max(Gabs(1:Nk,1:Nk),[],2);
            else
                [dddd,colMaxG_k]=max(Gabs,[],2);
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

function W=cca(X)
% W=cca(X)
%
% Canonical correlation analysis.
% Inputs:
%  X: Cell array of 2 entries (each entry is called a dataset).  Each cell entry contains
%  samples of dimension N1 x T and N2 x T.  Each dataset needs to have the same number of
%  samples, T.  But the dimensionality of each dataset can be different.
%

if nargin==0
    cca_test; % built-in test function
    return
end
X_iscell=true;
if ~iscell(X)
    X_iscell=false;
    Xmat=X;
    X=cell(2,1);
    X{1}=Xmat(:,:,1);
    X{2}=Xmat(:,:,2);
end

K=length(X);

if K~=2
    error('CCA only works for two datasets -- try calling multiset-CCA routines or an independent vector analysis routine.')
end
N1=size(X{1},1);
N2=size(X{2},1);
swapX=false;
if N2<N1
    swapX=true;
    temp=X{1};
    X{1}=X{2};
    X{2}=temp;
end
Rx=cell(K,K);
T=size(X{1},2);
for k1=1:K
    if size(X{k1},2)~=T
        error('Need same number of samples in each dataset.')
    end
    for k2=1:K
        Rx{k1,k2}=X{k1}*X{k2}'/T;
    end
end

temp=Rx{2,2}\Rx{2,1};
A=Rx{1,2}*temp;
[W1,dddd]=eig(A,Rx{1,1});
W2=temp*W1;
W=cell(K,1);
W{1}=vecnorm(W1)'; % transpose is to achieve our desired form
W{2}=vecnorm(W2)'; % vecnorm normalizes scaling of columns (arbitrary & unnecessary)

if swapX
    temp=W{1};
    W{1}=W{2};
    W{2}=temp;
end
if ~X_iscell
    Wcell=W;
    N=size(Wcell{1},1);
    W=zeros(N,N,K);
    W(:,:,1)=Wcell{1};
    W(:,:,2)=Wcell{2};
end
return


function Rx = computeCrossCov(X)

try
    
    C = size(X, 1);
    N = size(X, 3);
    X = permute(X, [2, 1, 3]);
    X = reshape(X, size(X, 1), size(X, 2)*size(X, 3));
    
    cov_m = X'*X/size(X, 1);
    
    Rx = mat2cell(cov_m, repmat(C, 1, N), repmat(C, 1, N));
    
catch
    
    [N, T, K] = size(X);
    
    Rx = cell(K,K);
    for k1=1:K
        Rx{k1,k1}=1/T*(X(:,:,k1)*X(:,:,k1)');
        for k2=(k1+1):K
            Rx{k1,k2}=1/T*(X(:,:,k1)*X(:,:,k2)');
            Rx{k2,k1}=Rx{k1,k2}';
        end
    end
    
end


function X = cholsol(S, B)
%% Cholesky decomposition

L = (chol(S, 'lower'));
Y = L \ B;
X = L' \ Y;