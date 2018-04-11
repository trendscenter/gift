function [Wwasobi,Wsobi,ISR,signals]= icatb_ewasobi(x,AR_order,rmax,eps0)
%
% implements algorithm WASOBI for blind source separation of
% AR sources in a fast way, allowing separation up to 100 sources
% in the running time of the order of tens of seconds.
%
% INPUT:    x .... input data matrix d x N
%           d .... dimension of the data
%           N .... length of the data
%           ARmax .. maximum AR order of the separated sources
%           rmax ... a constant that may help to stabilize the algorithm.
%                    it has the meaning of maximum magnitude of poles of the
%                    AR sources. The choice rmax=1 means that no stabilization
%                    is applied. The choice rmax=0.99 may lead to more stable
%                    results.
%           eps0 ... machine dependent constant to control condition number
%                    of weight matrices
%
% OUTPUT: Wwasobi ...... estimated de-mixing matrix
%         Wsobi ........ initial estimate of the matrix obtained by FFDiag2
%         ISR .......... estimated ISR matrix which represents approximate accuracy
%                        of the separation provided that there is no additive
%                        noise in the model.
%         signals....... separated signals
%
% Code by Petr Tichavsky, using inputs from Eran Doron, Pavel Laskov
% and Andreas Ziehe. A paper that describes the algorithm is submitted
% to ICASSP 2007.
%
if nargin<4
    eps0=5.0e-7;
end
if nargin<3
    rmax=0.99;
end
num_of_iterations = 3;
[d N]=size(x);
Xmean=mean(x,2);
x=x-Xmean*ones(1,N);  %%%%%%%%%  removing the sample mean
%[AOL_1 Rx_est] = sobi(x,AR_order); WOL_1=inv(AOL_1);
T=length(x(1,:))-AR_order;
C0=corr_est(x,T,AR_order);
%
Wsobi = ffdiag2(C0,eye(d));
Wwasobi=Wsobi; WOL_1=Wsobi;
Rx_est = [];
%
for in = 1:num_of_iterations
    [WOL_1,Rx_est,AR]=newasobi(x,AR_order+1,WOL_1,Rx_est,rmax,eps0);
    Wwasobi=WOL_1*Wwasobi;
end
ISR=CRLB4(AR)/N;
sqdnorm=diag(Wwasobi*C0(1:d,1:d)*Wwasobi');
Wwasobi=Wwasobi./(sqrt(sqdnorm)*ones(1,d));
signals=Wwasobi*x+(Wwasobi*Xmean)*ones(1,N);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of EWASOBI
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W_est_WASOBI,Rx_est,AR]=newasobi(x,M,W_pre,Rx_est1,rmax,eps0)
%
%   Computes one iteration needed for wasobi
%
%   d ... number of independent components
%   M ... number of lags = AR_order+1
%   Vnorm ... estimated demixing matrix obtained in preprocessing
%   Wnorm ...  estimated demixing matrix obtained by the advanced tuning
%

status=1;
[d,N]=size(x);
d2=d*(d-1)/2;
T=length(x(1,:))-M+1;
if isempty(Rx_est1)
    Rx_est=corr_est(x,T,M-1,W_pre);
else
    for index=1:M
        id=(index-1)*d;
        Rx_est(:,id+1:id+d)=W_pre*Rx_est1(:,id+1:id+d)*W_pre';
        %       Rx_est(:,:,index)=W_pre*Rx_est1(:,:,index)*W_pre';
    end
end
R=zeros(M,d);
for index=1:M
    id=(index-1)*d;
    R(index,:)=diag(Rx_est(:,id+1:id+d)).';
    %     R(index,:)=diag(Rx_est(:,:,index)).'; %%% columns of R will contain
    %%% covariance function of the separated components
end
%
[AR,sigmy]=armodel(R,rmax);      %%% compute AR models of estimated components
%
AR3=[];
for i=2:d
    for k=1:i-1
        AR3=[AR3 conv(AR(:,i),AR(:,k))];
    end
end
phi=ar2r(AR3);     %%%%%%%%%% functions phi to evaluate CVinv
CVinv=THinv5(phi,M,d2,eps0*phi(1,:));  %%%% to compute inversions of CV
%%%% It has dimension zeros(M,M*d2).
im=1;
for i=2:d
    for k=1:i-1
        fact=1/(sigmy(1,i)*sigmy(1,k));
        CVinv(:,(im-1)*M+1:im*M)=CVinv(:,(im-1)*M+1:im*M)*fact;
        im=im+1;
    end
end
W_est_WASOBI=serial10(Rx_est, CVinv);  %%% TO MINIMIZE THE WEIGHTED CRITERION

%%%%%%%%%%%%%%%%%%%%%%%%%%   of newasobi
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ r ] = ar2r( a )
%%%%%
%%%%% Computes covariance function of AR processes from
%%%%% the autoregressive coefficients using an inverse Schur algorithm
%%%%% and an inverse Levinson algorithm (for one column it is equivalent to
%%%%%      "rlevinson.m" in matlab)
%
if (size(a,1)==1)
    a=a'; % chci to jako sloupce
end

[p m] = size(a);    % pocet vektoru koef.AR modelu
alfa = a;
K=zeros(p,m);
p = p-1;
for n=p:-1:1
    K(n,:) = -a(n+1,:);
    for k=1:n-1
        alfa(k+1,:) = (a(k+1,:)+K(n,:).*a(n-k+1,:))./(1-K(n,:).^2);
    end
    a=alfa;
end
%
r  = 1./prod(1-K.^2);
f=[r(1,:); zeros(p,m)];
b=f;
for k=1:p
    for n=k:-1:1
        f(n,:)=f(n+1,:)+K(n,:).*b(k-n+1,:);
        b(k-n+1,:)=-K(n,:).*f(n+1,:)+(1-K(n,:).^2).*b(k-n+1,:);
    end
    b(k+1,:)=f(1,:);
    r=[r; f(1,:)];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%  of ar2r
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function R_est=corr_est(x,T,q,W_pre)
%
[K,KK]=size(x);
R_est=[];
if ~exist('W_pre')
    for index=1:q+1
        R_est=[R_est 1/T*(x(:,1:T)*x(:,index:T+index-1)')];
    end
else
    for index=1:q+1
        R_est=[R_est 1/T*W_pre*(x(:,1:T)*x(:,index:T+index-1)')*W_pre'];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of corr_est
%
%xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx


function [AR,sigmy]=armodel(R,rmax)
%
% to compute AR coefficients of the sources given covarience functions
% but if the zeros have magnitude > rmax, the zeros are pushed back.
%
[M,d]=size(R); v1=zeros(1,d); v2=v1; Rs=R;
for id=1:d
    AR(:,id)=[1; -toeplitz(R(1:M-1,id),R(1:M-1,id)')\R(2:M,id)];
    v=roots(AR(:,id)); %%% mimicks the matlab function "polystab"
    %    v1(1,id)=max(abs(v));
    vs=0.5*(sign(abs(v)-1)+1);
    v=(1-vs).*v+vs./conj(v);
    vmax=max(abs(v));
    %    v2(1,id)=max(abs(v));
    if vmax>rmax
        v=v*rmax/vmax;
    end
    AR(:,id)=real(poly(v)'); %%% reconstructs back the covariance function
end
Rs=ar2r(AR);
sigmy=R(1,:)./Rs(1,:);
% [v1; v2]
%%%%%%%%%%%%%%%%%%%%%%%  of armodel

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function W=serial10(Rx_est,CVinv)
%
[d Md]=size(Rx_est);
M=Md/d; dd=d^2-d; dd2=dd/2;
W=eye(d);
%
num_of_iteration=5;
for iteration=1:num_of_iteration
    A0=eye(d);
    im=1;
    for i=2:d
        for k=1:i-1
            Y(im,:)=(Rx_est(i,k:d:Md)+Rx_est(k,i:d:Md))/2;
            %        Y(im,:)=squeeze(Rx_est(i,k,:)+Rx_est(k,i,:))/2;
            im=im+1;
        end
    end
    for i=1:d
        LaMat(i,:)=Rx_est(i,i:d:Md);
    end
    b11=zeros(dd2,1); b12=b11; b22=b11; c1=b11; c2=c1;
    m=0;
    for id=2:d
        for id2=1:id-1
            m=m+1; im=(m-1)*M;
            Wm=CVinv(:,im+1:im+M);
            Wlam1=Wm*LaMat(id,:)';
            Wlam2=Wm*LaMat(id2,:)';
            b11(m,1)=LaMat(id2,:)*Wlam2;
            b12(m,1)=LaMat(id,:)*Wlam2;
            b22(m,1)=LaMat(id,:)*Wlam1;
            c1(m,1)=c1(m,1)+Wlam2'*Y(m,:)';
            c2(m,1)=c2(m,1)+Wlam1'*Y(m,:)';
        end
    end
    det0=b11.*b22-b12.^2;
    d1=(c1.*b22-b12.*c2)./det0;
    d2=(b11.*c2-b12.*c1)./det0;
    %    value=norm([d1; d2])
    m=0;
    for id=2:d
        for id2=1:id-1
            m=m+1;
            A0(id,id2)=d1(m,1);
            A0(id2,id)=d2(m,1);
        end
    end
    Ainv=inv(A0);
    if iteration<num_of_iteration
        for m=1:M
            %          Rx_est(:,:,m)=Ainv*Rx_est(:,:,m)*Ainv';
            id=(m-1)*d;
            Rx_est(:,id+1:id+d)=Ainv*Rx_est(:,id+1:id+d)*Ainv';
        end
    end
    W=Ainv*W;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  of serial10


function ISR = CRLB4(AR)
%
% CRLB4(AR) generates the CRLB for gain matrix elements (in term
% of ISR) for blind separation of K Gaussian autoregressive sources
% whose AR coefficients (of the length M, where M-1 is the AR order)
% are stored as columns in matrix AR.

[M K]=size(AR);

Rs=ar2r(AR);

sum_Rs_s=zeros(K,K);

for s=0:M-1
    for t=0:M-1
        sum_Rs_s=sum_Rs_s+(AR(s+1,:).*AR(t+1,:))'*Rs(abs(s-t)+1,:);
    end
end

denom=sum_Rs_s'.*sum_Rs_s+eye(K)-1;
ISR=sum_Rs_s'./denom.*(ones(K,1)*Rs(1,:))./(Rs(1,:)'*ones(1,K));
ISR(eye(K)==1)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of CRLB4

function G=THinv5(phi,K,M,eps)
%
%%%% Implements fast (complexity O(M*K^2))
%%%% computation of the following piece of code:
%
%C=[];
%for im=1:M
%  A=toeplitz(phi(1:K,im),phi(1:K,im)')+hankel(phi(1:K,im),phi(K:2*K-1,im)')+eps(im)*eye(K);
%  C=[C inv(A)];
%end
%
% DEFAULT PARAMETERS: M=2; phi=randn(2*K-1,M); eps=randn(1,2);
%   SIZE of phi SHOULD BE (2*K-1,M).
%   SIZE of eps SHOULD BE (1,M).

phi(2*K,1:M)=0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
almold=2*phi(1,:)+eps;
C0=1./almold;
x1=zeros(K,M); x2=x1; x3=x1; x4=x1;
x1(1,:)=C0; x2(1,:)=C0;
x3(1,:)=-C0.*phi(2,:);
x4(1,:)=-2*C0.*phi(2,:);
x4old=[];
lalold=2*phi(2,:)./almold;
for k=1:K-1
    f2o=phi(k+1:-1:2,:)+phi(k+1:2*k,:);
    alm=sum(f2o.*x4(1:k,:),1)+phi(1,:)+eps+phi(2*k+1,:);
    a0=zeros(1,M);
    if k<K-1
        a0=phi(k+2,:);
    end
    gam1=sum(f2o.*x1(1:k,:),1);
    gam3=sum(f2o.*x3(1:k,:),1)+a0+phi(k,:);
    x4(k+1,:)=ones(1,M);
    b1m=sum(([phi(2:k+1,:); a0]+[zeros(1,M); phi(1:k,:)]).*x4(1:k+1,:));
    b2m=sum(([a0; phi(k+1:-1:2,:)]+phi(k+2:2*k+2,:)).*x4(1:k+1,:));
    latemp=b2m./alm;
    b2m=latemp-lalold; lalold=latemp;
    bom=alm./almold;
    ok=ones(k+1,1);
    x2(1:k+1,:)=x4(1:k+1,:).*(ok*(1./alm));
    x1(1:k+1,:)=[x1(1:k,:); zeros(1,M)]-(ok*gam1).*x2(1:k+1,:);
    x3(1:k+1,:)=[x3(1:k,:); zeros(1,M)]-(ok*gam3).*x2(1:k+1,:);
    x4temp=x4(1:k,:);
    x4(1:k+1,:)=[zeros(1,M); x4(1:k,:)]+[x4(2:k,:); ones(1,M); zeros(1,M)]...
        -(ok*bom).*[x4old; ones(1,M); zeros(1,M)]...
        -(ok*b2m).*x4(1:k+1,:)-(ok*b1m).*x1(1:k+1,:)-(ok*x4(1,:)).*x3(1:k+1,:);
    x4old=x4temp;
    almold=alm;
end
MK=M*K;
G=zeros(K,MK);
G(:,1:K:MK)=x1; clast=zeros(K,M);
f1=[phi(2:K,:); zeros(1,M)]+[zeros(1,M); phi(1:K-1,:)];
f2=[zeros(1,M); phi(K:-1:2,:)]+[phi(K+1:2*K-1,:); zeros(1,M)];
for k=2:K
    ck=G(:,k-1:K:MK);
    G(:,k:K:MK)=[ck(2:K,:); zeros(1,M)]+[zeros(1,M);  ck(1:K-1,:)]...
        -clast-(ok*sum(f1.*ck)).*x1-(ok*sum(f2.*ck)).*x2-(ok*ck(1,:)).*x3...
        -(ok*ck(K,:)).*x4;
    clast=ck;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   of THinv5

function [W] = ffdiag2(Rx_est,ini)
%
% computes Ziehe's ffdiag with pre-processing for maintaining
% high stability even for badly conditioned mixing matrices
% Next, ffdiag requires matrices stored in 3-dimensional field.
[d Md]= size(Rx_est);
M=Md/d;
[H1,D] = eig(Rx_est(:,1:d));
W_est = [];
for k = 1:d
    W_est = [W_est; D(k,k)^(-.5)*H1(:,k)'/(H1(:,k)'*H1(:,k))];
end
R_zest=zeros(d,d,M);
for k = 1:M
    id=(k-1)*d;
    R_zest(:,:,k) = W_est*Rx_est(:,id+1:id+d)*W_est';
end
%%%
[W C stat] = ffdiag(R_zest,ini);  %%% CALLING ORIGINAL FFDIAG
%%%
W=W*W_est;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   of ffdiag2

function [V,C,stat] = ffdiag(C0,V0)
%FFDIAG  Diagonalizes a set of matrices C_k with a (single) non-orthogonal transformation.
%
%
% Usage:
%
%  function [V,CD,stat] = ffdiag(C0,V0)
%
% Input:
%    C0 - N x N x K array containing K symmetric matrices C_k which are to be diagonalized.
%    V0 - initial value for the diagonalizer (default eye(N))
%
% Output:
%    V    - joint diagonalization transformation (dimension N x N)
%    CD   - diagonalized set of matrices (dimension N x N x K).
%    stat - structure with some statistics of the algorithm's performance
%        .etime      - elapsed time
%        .errdiag    - diagonalization error
%
% code by Andreas Ziehe and Pavel Laskov
%
%  (c) 2004  Fraunhofer FIRST.IDA
%
% Usage example:
%
% gen_mat; [V,C,stat]=ffdiag(C0); imagesc(V*A);



%
% The algorithm is derived in the paper:
%  A. Ziehe, P. Laskov, G. Nolte and K.-R. Mueller,
% "A Fast Algorithm for Joint Diagonalization with Non-orthogonal
%  Transformations and its Application to Blind Source Separation"
%  Journal of Machine Learning Research vol 5, pages 777-800, 2004.
%
% An earlier version has been presented at the ICA 2003 Workshop in
% Nara, Japan
%
%  A. Ziehe and P. Laskov and K.-R. Mueller and G. Nolte
%  "A Linear Least-Squares Algorithm for Joint Diagonalization",
%  in Proc. of the 4th International Symposium on
%  Independent Component Analysis and Blind Signal Separation
%  (ICA2003), pages 469--474, Nara, Japan 2003.
%
%
% See also
% http://www.first.fraunhofer.de/~ziehe/research/ffdiag.html


eps = 1e-9;

%theta = 0.1;    % threshold for stepsize heuristics

[m,n,k] = size(C0);
C = C0;
Id=eye(n);
V = Id;

%more clever initialization ???
%[V,D]=eig(C0(:,:,1),C0(:,:,2));
%or
%[V,D]=eig(sum(C0,3));
%V=V';

inum = 1;
df = 1;

stat.method='FFDIAG';
tic;
while (df > eps & inum < 100)

    % Compute W
    W = getW(C);

    %old normalization  no longer recommended
    %as it deteriorates  convergence
    %
    %if (norm(W,'fro') > theta)
    %  W = theta/norm(W,'fro')*W;
    %end

    % A much better way is to
    % scale W by power of 2 so that its norm is <1 .
    % necessary to make approximation assumptions hold
    % i.e. W should be small
    %

    %if 0,
    [f,e] = log2(norm(W,'inf'));
    % s = max(0,e/2);

    s = max(0,e-1);
    W = W/(2^s );

    % Compute update
    V = (Id+W)*V;

    %re-normalization
    V=diag(1./sqrt(diag(V*V')))*V;  %norm(V)=1

    for i=1:k
        C(:,:,i) = V*C0(:,:,i)*V';
    end

    % Save stats
    stat.W(:,:,inum) = W;
    [f] = get_off(V,C0);

    stat.f(inum) = f;
    stat.errdiag(inum)=cost_off(C0,normit(V')');
    stat.nw(inum) = norm(W(:));


    if inum > 2
        df = abs(stat.f(inum-1)-stat.f(inum));

    end

    %  fprintf(1,'itr %d :  %6.8f  \n ',inum, f);
    inum = inum+1;
end

stat.etime=toc;   %elapsed time
stat.niter=inum;  %number of iterations
stat.V=V;         %estimated diagonalizer

%subplot(2,1,1);
%semilogy(stat.f);
%title('objective function');

%subplot(2,1,2);
%semilogy(stat.nw');
%title('norm W');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of ffdiag

function W = getW(C)
%getW   computes the update for FFDIAG.
%
% Usage:
%    function W = getW(C)
%
% Input:
%    C -  3-d array of matrices C_k
%
%  Output:
%    W - update matrix with zeros on the main diagonal
%
% See also:
%          FFDIAG
%


[m,n,K] = size(C);
W = zeros(n);
z = zeros(n);
y = zeros(n);

for i=1:n
    for j=1:n

        for k=1:K
            z(i,j) = z(i,j)+C(i,i,k)*C(j,j,k);
            y(i,j) = y(i,j)+0.5*C(j,j,k)*(C(i,j,k)+conj(C(j,i,k)));
        end
    end
end

for i=1:n-1
    for j=i+1:n
        W(i,j) = (z(j,i)*y(j,i)-z(i,i)*y(i,j))/(z(j,j)*z(i,i)-z(i,j)^2);
        W(j,i) = ((z(i,j)*y(i,j)-z(j,j)*y(j,i))/(z(j,j)*z(i,i)-z(i,j)^2));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of getW

function V=normit(V)
%NORMIT  Normalize the Frobenius norm of all columns of a matrix to 1.
[n,m]=size(V);

for k=1:n,
    nn=norm(V(:,k),'fro');
    if nn<eps,warning('division may cause numerical errors');end
    V(:,k)=V(:,k)/nn;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of normit

function [cost]=cost_off(C,V)
%COST_OFF computes diagonalization error as the norm of the off-diagonal elements.
%  C set of matrices in a 3-d array
%  V diagonalizing matrix


[N,N,K]=size(C);

if nargin>1,
    for k=1:K,
        C(:,:,k)=V*C(:,:,k)*V';
    end
end

cost=0;
for k=1:K
    Ck=C(:,:,k);Ck=Ck(:);
    Ck(1:N+1:N*N)=0;
    cost=cost+norm(Ck)^2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of cost_off

function [f,g] = get_off(V,C)
%GET_OFF  a wrapper for OFF
%
% This function loops over 'off' and  adds up things over multiple matrices C.

[m,n,K] = size(C);
f = 0;
g = sparse(n,n);
h = sparse(n^2,n^2);
h1 = sparse(n^2,n^2);
h2 = sparse(n^2,n^2);

for k=1:K
    if nargout <= 1
        f = f + off(V,C(:,:,k));
    elseif nargout == 2
        [ff,gg] = off(V,C(:,:,k));
        f = f+ff;
        g = g+gg;
        gg
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of get_off

function [f,g] = off(V,C)
%OFF - computes the value and the gradient of the "off" diagonality measure
%
% Usage:
%        [f,g] = off(V,C)
%
% code by Pavel Laskov and Andreas Ziehe
%
% (c) 2004 Fraunfofer FIRST.IDA

F = V*C*V';
f = trace(F'*F) - trace(F.*F);

if nargout > 1
    g = 4*(F-diag(diag(F)))*V*C;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% of off  %%% and of the whole ffdiag