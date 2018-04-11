function [out_M, W] = icatb_multi_fixed_ICA_R_Cor(mixedsig,refsig,loopnum, threshold,i_Rou,i_Uk,i_Gama,epsilon)
%Fixed point multi-unit ICA-R with Cor measure
%
% update: 2006.9.6

% mixedsig   the mixed signal
% refsig     reference signal
% loopnum    loop number for w learning
% threshold  for the constraint g(w)
% i_Rou      positive constant
% i_Uk       Lagrange multiplier
% i_Gama     penelty parameter
% epsilon    learning error

% out_M      extracted signal
% W          optimum W

tic;
% Removes the mean of mixedsig and refsig.
[mixedsigr,meanmixed]=icatb_remmean(mixedsig);
refsigr=icatb_remmean(refsig);
%------- End of mean removing ---------% 

%Pre-whitening mixedsig and refsig
Rxx = cov(mixedsigr');
[AV,B] = eig(Rxx);
PreWhiten = inv(sqrt(B))*AV';
DepreWhiten=AV*sqrt(B);
X = PreWhiten*mixedsigr;
%
r = diag(1./std(refsigr,0,2))*refsigr;
%-------End of Pre-whitening------%

%_____________Parameters Setting______________%
[dim,numsamp] = size(mixedsig); % number of sample
[dimr,numsampr]=size(r);

threshold=threshold*ones(dimr,1);

Rou= i_Rou*ones(dimr,1);    % positive constant in the contrast function, 
Uk=i_Uk*ones(dimr,1);       % Lagrangian multiplier for g(w), recording all Uk learned,
Gama = i_Gama*ones(dimr,1); % scalar penalty parameter; 

v=randn(dimr,numsamp);     % Ganssian variable with zero mean and unity variance
W=zeros(dim,dimr);         % weight vector, recording all w learned
Wold=W; 
y=zeros(dimr,numsamp);     % output signal, recording all y estimated

for j=1:1:loopnum;
    %
    W=X*(icatb_Gdfun(W'*X))'*diag(Rou)/numsamp+0.5*X*r'*diag(Uk)/numsamp; %%%%%Cor
    W= W/norm(W);% normalizing W
    y=W'*X;
    %   
    for k=1:dimr; %%%%%Cor
        if Uk(k)+Gama(k)*icatb_gw_yr_Cor(W(:,k)'* X,r(k,:),numsamp,threshold(k))>0;
           Uk(k)=Uk(k)+Gama(k)*icatb_gw_yr_Cor(W(:,k)'* X,r(k,:),numsamp,threshold(k));
        else
           Uk(k)=0;
        end
    end    
    
    Rou = mean(icatb_Gfun(y),2)-mean(icatb_Gfun(v),2);
    
    %---------Decorrelation of W---------
    W=W';
    W1=sqrtm(W*W');
    W2=inv(W1);
    W=W2*W;
    W=W';
    %---------end of decorrelation of W---------
    
    minAbsW = min(abs(diag(W' * Wold)));
    if (1 - minAbsW < epsilon)
           fprintf('computed ( %d steps ) \n', j);
           break;
    end
    Wold=W;
end
t_fixed_Cor=toc
%%%%%%-----end of Fixed point ICA-R------%%%%%%

% % outputs
W=W'*PreWhiten;
out_M=y;
