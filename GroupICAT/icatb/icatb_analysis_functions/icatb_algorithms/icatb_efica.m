function [Wefica, ISRef, Hef, Wsymm, ISRsymm, status]= icatb_efica(X, ini, FASTICApars, Saddlepars,EFICApars)
% function [Wefica, ISRef, Hef, Wsymm, ISRsymm, status]=efica(X, ini, g, SaddleTest)

%EFICA: [Wefica, ISRef, Hef, Wsymm, SIRsymm, status]=efica(X, ini, g, SaddleTest)
%
% version: 1.9  release: 24.11.2006
%
%          %improved speed
%
%Input data:
% X ... mixed data dim x N, where dim is number of signals and N is number of
%       samples
% ini ... starting point of the iterations for W matrix
% g ... nonlinearity used in the symmetric algorithm ('tanh'-'rati'-'gaus'-'pow3')
% SaddleTest ... if true (default) the test of saddle points on the demixed signals is done
%
%
%Output data:
% W_efica - demixing matrix produced by EFICA
% ISRef - ISR matrix estimator of EFICA components
% Hef - matrix for EFICA bias computation (in the noisy environment)
% Wsymm - demixing matrix produced by symmetric Fast-ICA
% ISRsymm - ISR matrix estimator of Fast-ICA components
% status - one if test of saddle points was positive, zero otherwise
%
% Examples of usage:
%
%  [Wefica, ISRef, Hef, Wsymm, ISRsymm]=efica(X);
%
%    use default settings (random initialization, 'tanh' nonlinearity,
%    the test of saddle points will be done)
%
%
%  [Wefica, ISRef, Hef, Wsymm, SIRsymm]=efica(X, randn(dim), 'rati', true);
%
%    random initialization, 'rati' nonlinearity used in symmetric FastICA,
%    the test of saddle points will be done
%
%
%  [Wefica, ISRef, Hef, Wsymm, SIRsymm]=efica(X, eye(dim), 'tanh', false);
%
%    random initialization, 'tanh' nonlinearity used in symmetric FastICA,
%    without the test of saddle points
%
%
%
%
%

[dim N]=size(X);

%Default values of parameters
if nargin<5 %EFICApars
    %     EFICApars = true;
    EFICAMaxIt=50; %%Maximum number of improving iterations
    eficaGaussianNL=3; %'ggda';  %Nonlinearity used for eficagaussian signals 'gaus'-'ggda'-'npnl'
else
    EFICAMaxIt=EFICApars.EFICAMaxIt ; %%Maximum number of improving iterations
    eficaGaussianNL=EFICApars.eficaGaussianNL;  %Nonlinearity used for eficagaussian signals 'gaus'-'ggda'-'npnl'
end

if nargin<4     %Saddlepars
    SaddleTest = 1;
    test_of_saddle_points_nonln = 4 ; %'rtsp';
    MaxItAfterSaddleTest=30; %%Maximum number of iterations after a saddle point was indicated
else
    SaddleTest = Saddlepars.SaddleTest ;
    test_of_saddle_points_nonln=Saddlepars.test_of_saddle_points_nonln;
    MaxItAfterSaddleTest=Saddlepars.MaxItAfterSaddleTest; %%Maximum number of iterations after a saddle point was indicated

end
if nargin<3 %FASTICApars
    g= 3; %'rati';
    MaxIt=100; %% Maximum number of FastICA iterations
else
    g=FASTICApars.g;
    MaxIt=FASTICApars.MaxIt; %% Maximum number of FastICA iterations

end

if nargin<2

    while 1
        ini=randn(dim);
        if cond(ini) < 10^3
            break
        end
    end
end

epsilon=0.0001; %Stop criterion
fineepsilon=1e-5; %Stop criterion for post-estimation
repeat=1;
rot2d=[1/sqrt(2) 1/sqrt(2);-1/sqrt(2) 1/sqrt(2)];
status=0;
min_correlation=0.8; %additive noise...0.75, noise free... 0.95, turn off (unstable)...0

%removing mean
X = X-mean(X,2)*ones(1,N);


%fprintf('EFICA...\n');

%preprocessing
C = cov(X');
CC = C^(-1/2);
Z = CC*X;

%FastICA
W=ini;
W_before_decor=W;
W=W*real(inv(W'*W)^(1/2));
%Wold=zeros(dim,dim);
crit=zeros(1,dim);
NumIt=0;
while repeat
    %%% Symmetric approach
    while (1-min(crit)>epsilon && NumIt<MaxIt)% && sum(double(changed>10))<2)
        Wold=W;
        switch g
            case 1 %'tanh'
                hypTan = tanh(Z'*W);
                W=Z*hypTan/N-ones(dim,1)*sum(1-hypTan.^2).*W/N;
            case 2 %'pow3'
                W=(Z*((Z'*W).^ 3))/N-3*W;
            case 3 %'rati'
                U=Z'*W;
                Usquared=U.^2;
                RR=4./(4+Usquared);
                Rati=U.*RR;
                Rati2=Rati.^2;
                dRati=RR-Rati2/2;
                nu=mean(dRati);
                %    beta=mean(Rati2);
                hlp=Z*Rati/N;
                %     mu=diag(W'*hlp)';

                %      J=ones(1,dim); tau=abs(nu-mu); gam=(beta-mu.^2);
                %      ISR=(gam'*J+J'*gam+J'*tau.^2)./(tau'*J+J'*tau).^2;
                %      ISR=ISR-diag(diag(ISR));

                W=hlp-ones(dim,1)*nu.*W;

            case 4%'gaus'
                U=Z'*W;
                Usquared=U.^2;
                ex=exp(-Usquared/2);
                gauss=U.*ex;
                dGauss=(1-Usquared).*ex;
                W=Z*gauss/N-ones(dim,1)*sum(dGauss).*W/N;
        end
        %   crit=abs(sum((W*diag(sqrt(1./sum(W.^2)))).*(W_before_decor*diag(sqrt(1./sum(W_before_decor.^2))))));
        W_before_decor=W;
        W=W*real(inv(W'*W)^(1/2));
        crit=abs(sum(W.*Wold));
        NumIt=NumIt+1;
    end %while iteration
    fprintf('EFICA v1.9 NumIt: %d\n',NumIt);
    repeat=0;
    %%%The SaddleTest of the separated components
    if SaddleTest
        SaddleTest=false; %%The SaddleTest may be done only one times
        u=Z'*W;
        switch test_of_saddle_points_nonln
            case 1 %'tanh'
                table1=(mean(log(cosh(u)))-0.37456).^2;
            case 2  %'gaus'
                table1=(mean(u.*exp(-u.^2))-1/sqrt(2)).^2;
            case 3 %'rati'
                table1=(mean(2*log(4+u.^2))-3.16).^2;
            case 4 %'rtsp'
                table1=(mean(u.^2./(1+abs(u)))-0.4125).^2;
            case 5 %'pow3'
                table1=(mean((pwr(u,4)))-3).^2;
        end
        rotated=zeros(1,dim);
        checked=1:dim;
        for i=checked
            for j=checked(checked>i)
                if (~rotated(i) && ~rotated(j))
                    h=[u(:,i) u(:,j)]*rot2d;   % rotate two columns i,j
                    switch test_of_saddle_points_nonln
                        case 1 %'tanh'
                            ctrl=(mean(log(cosh(h)))-0.37456).^2;
                        case 2 %'gaus'
                            ctrl=(mean(exp(-h.^2/2)-1/sqrt(2))).^2;
                        case 3 %'rati'
                            ctrl=(mean(2*log(4+h.^2))-3.16).^2;
                        case 4 %'rtsp'
                            ctrl=(mean(h.^2./(1+abs(h)))-0.4125).^2;
                        case 5 %'pow3'
                            ctrl=(mean((pwr(h,4)))-3).^2;
                    end
                    if sum(ctrl)>table1(i)+table1(j)
                        %bad extrem indicated
                        rotated(i)=1;rotated(j)=1; %do not test the rotated signals anymore
                        W(:,[i j])=W(:,[i j])*rot2d;
                        repeat=1; %continue in iterating - the test of saddle points is positive
                        NumIt=MaxIt-MaxItAfterSaddleTest;
                        fprintf('EFICA: rotating components: %d and %d\n',i,j);
                        status=1;
                    end
                end
            end
        end
    end %if SaddleTest
    crit=zeros(1,dim);
end %while repeat

Wsymm=W'*CC;

%estimated signals
s=W'*Z;

%estimate SIRs of the symmetric approach
switch g
    case  1 %'tanh'
        mu=mean(s.*tanh(s),2);
        nu=mean(1./cosh(s).^2,2);
        beta=mean(tanh(s).^2,2);
    case 2 %'rati'
        ssquare=s.^2;
        mu=mean(ssquare./(1+ssquare/4),2);
        nu=mean((1-ssquare/4)./(ssquare/4+1).^2,2);
        beta=mean(ssquare./(1+ssquare/4).^2,2);
    case 3 %'gaus'
        aexp=exp(-s.^2/2);
        mu=mean(s.^2.*aexp,2);
        nu=mean((1-s.^2).*aexp,2);
        beta=mean((s.*aexp).^2,2);
    case 4 %'pow3'
        mu=mean(s.^4,2);
        nu=3*ones(dim,1);
        beta=mean(s.^6,2);
end
J=ones(1,dim); gam=(nu-mu)'; jm=(beta-mu.^2)';
Werr=(jm'*J+J'*jm+J'*gam.^2)./(abs(gam)'*J+J'*abs(gam)).^2;
Werr=Werr-diag(diag(Werr));
ISRsymm=Werr/N;
%SIRsymm=-10*log10(sum(Werr,2)'/N);

ekurt=mean(pwr(s,4),2); %%% estimate fourth moment

%Werr=Werr/N+diag((ekurt-1)/(4*N));
% GCov=diag(Werr(:));
% for i=1:dim-1
%     for j=i+1:dim
%         Gjiij=(-jm(i)-jm(j)+abs(gam(i)*gam(j)))/(abs(gam(i))+abs(gam(j)))^2;
%         GCov((i-1)*dim+j,(j-1)*dim+i)=Gjiij;
%         GCov((j-1)*dim+i,(i-1)*dim+j)=Gjiij;
%     end
% end
%varW_symm=reshape(diag(kron(Wsymm',eye(dim))*GCov*kron(Wsymm,eye(dim))),dim,dim);

%varA_symm=reshape(diag(kron(eye(dim),inv(Wsymm))*GCov*kron(eye(dim),inv(Wsymm)')),dim,dim);

%EFICA


% if eficaGaussianNL=='npnl'
%     ekurt=ones(1,dim)*4;
% end
for j=1:dim
    w=W(:,j);
    if ekurt(j)<2.4184   %%% sub-Gaussian -> try "GGD score function" alpha=>3
        if ekurt(j)<1.7  %%% the distribution seems to be extremly subgaussian
            %alpha=50; % bimodal-gaussian distribution will be considered
            alpha=15;
        else %%% estimate shape parameter alpha from the fourth moment
            if ekurt(j)>1.8
                alpha=1/(sqrt((5*ekurt(j)-9)/6/pi^2)+(1.202)*3*(5*ekurt(j)-9)/pi^4);
            else
                alpha=15; %the distribution is likely uniform
            end
        end
        if alpha<3 %%% the distribution is rather gaussian -> keep the original result
            w=W_before_decor(:,j);
        elseif alpha<=15 %%% try score function sign(x)*|x|^(alpha-1)
            wold=zeros(dim,1);
            nit=0;
            alpha=ceil(min([alpha 15]));
            while ((1-abs(w'*wold/norm(w))>fineepsilon) && (nit<EFICAMaxIt) &&...
                    (abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))>min_correlation))
                w=w/norm(w);
                %abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))
                wold=w;
                u=Z'*w;
                ualpha=pwr(abs(u),alpha-2);
                w=Z*(u.*ualpha)/N-(alpha-1)*mean(ualpha)*w;
                nit=nit+1;
            end
            sest=(w/norm(w))'*Z;
            if abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))>min_correlation,
                mu(j)=mean(pwr(abs(sest),alpha));
                nu(j)=(alpha-1)*mean(pwr(abs(sest),alpha-2));
                beta(j)=mean(pwr(abs(sest),2*alpha-2));
            end
        else  %trying bimodal distribution
            wold=zeros(dim,1);
            nit=0;
            while ((1-abs(w'*wold/norm(w))>fineepsilon) && (nit<EFICAMaxIt) &&...
                    (abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))>min_correlation))
                w=w/norm(w);
                wold=w;
                u=Z'*w;
                m=mean(abs(u)); %estimate the distance of distribution's maximas
                e=sqrt(1-m^2); %then their variance is..., because the overall variance is 1
                if e<=0.05, e=0.05; m=sqrt(1-e^2); end %due to stability
                uplus=u+m; uminus=u-m;
                expplus=exp(-uplus.^2/2/e^2);
                expminus=exp(-uminus.^2/2/e^2);
                expb=exp(-(u.^2+m^2)/e^2);
                g=-(uminus.*expminus + uplus.*expplus)./(expplus+expminus)/e^2;
                gprime=-(e^2*(expplus.^2+expminus.^2)+(2*e^2-4*m^2)*expb)./(expplus+expminus).^2/e^4;
                w=Z*g/N-mean(gprime)*w;
                nit=nit+1;
            end
            u=(w/norm(w))'*Z;
            if abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))>min_correlation,
                m=mean(abs(u)); e=sqrt(1-m^2);
                uplus=u+m; uminus=u-m;
                expplus=exp(-uplus.^2/2/e^2); expminus=exp(-uminus.^2/2/e^2); expb=exp(-(u.^2+m^2)/e^2);
                g=-(uminus.*expminus + uplus.*expplus)./(expplus+expminus)/e^2;
                gprime=-(e^2*(expplus.^2+expminus.^2)+(2*e^2-4*m^2)*expb)./(expplus+expminus).^2/e^4;
                mu(j)=mean(u.*g);
                nu(j)=mean(gprime);
                beta(j)=mean(g.^2);
            end
        end %bi-modal variant
    elseif ekurt(j)>3.762      %%% supergaussian (alpha<1.5)
        switch eficaGaussianNL
            case 1 %'npnl' %nonparametric density model; the same as in NPICA algorithm; very slow
                wold=zeros(dim,1);
                nit=0;
                while ((1-abs(w'*wold/norm(w))>fineepsilon) && (nit<EFICAMaxIt) &&...
                        (abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))>min_correlation))
                    w=w/norm(w);
                    [G GP]=nonln(w'*Z);
                    wold=w;
                    w=Z*G'/N-mean(GP)*w;
                    nit=nit+1;
                end
                [G GP]=nonln((w/norm(w))'*Z);
                mu(j)=mean(((w/norm(w))'*Z).*G);nu(j)=mean(GP);beta(j)=mean(G.^2);
            case 2 %'gaus'
                wold=zeros(dim,1);
                nit=0;
                while ((1-abs(w'*wold/norm(w))>fineepsilon) && (nit<EFICAMaxIt) &&...
                        (abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))>min_correlation))
                    w=w/norm(w);
                    aexp=exp(-(w'*Z).^2/2);
                    wold=w;
                    w=Z*((w'*Z).*aexp)'/N-w*mean((1-(w'*Z).^2).*aexp);
                    nit=nit+1;
                end
                sest=(w/norm(w))'*Z; aexp=exp(-sest.^2/2);
                mu(j)=mean(sest.^2.*aexp); nu(j)=mean((1-sest.^2).*aexp);
                beta(j)=mean((sest.*aexp).^2);
            case 3 %'ggda' %our proposal
                gam=3.3476;
                wold=zeros(dim,1);
                nit=0;
                while ((1-abs(w'*wold/norm(w))>fineepsilon) && (nit<EFICAMaxIt) &&...
                        (abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))>min_correlation))
                    w=w/norm(w);
                    u=w'*Z;
                    ualpha=abs(u);
                    aexp=exp(-gam*ualpha);
                    wold=w;
                    w=Z*(u.*aexp)'/N-w*mean((1-gam*ualpha).*aexp);
                    nit=nit+1;
                end
                if abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))>min_correlation,
                    sest=(w/norm(w))'*Z; ualpha=abs(sest); aexp=exp(-gam*ualpha);
                    mu(j)=mean(sest.^2.*aexp); nu(j)=mean((1-gam*ualpha).*aexp);
                    beta(j)=mean((sest.*aexp).^2);
                end
        end
    else  % alpha \in (1.5, 3) ; keep the original nonlinearity
        w=W_before_decor(:,j);
    end
    if abs((W(:,j)/norm(W(:,j)))'*(w/norm(w)))<min_correlation
        %%%%  the signal changed too much, thus,
        %%%%  seems to be rather gaussian -> keep the original
        %%%%  nonlinearity and result
        W(:,j)=W_before_decor(:,j);
    else
        W(:,j)=w;
    end
end

%Refinement
J=abs(mu-nu)';
B=(beta-mu.^2)';
%EJ=abs(Smu-Snu)';
%EB=(Sbeta-Smu.^2)';
Werr=zeros(dim); Hef=zeros(dim);
for k=1:dim
    ccc=(J*B(k)/J(k))./(B+J.^2);ccc(k)=1;
    Hef(k,:)=(ccc.*J-ccc(k)*J(k))./(ccc.*J+ccc(k)*J(k));
    WW=W*diag(ccc);
    sloupce=setdiff(1:dim,find(sum(WW.^2)/max(sum(WW.^2))<1e-7)); % removing almost zero rows
    M=WW(:,sloupce);
    M=M*real(inv(M'*M)^(1/2));
    if sum(sloupce==k)==1
        Wefica(:,k)=M(:,sloupce==k);
    else
        w=null(M');
        if size(w,2)==1
            Wefica(:,k)=w(:,1);
        else % there are probably more than two gaussian components => use the old result to get a regular matrix
            Wefica(:,k)=Wsymm(:,k);
        end
    end
    %Estimate variance of elements of the gain matrix
    Werr(k,:)=B(k)*(B+J.^2)./(J.^2*B(k)+J(k)^2*(B+J.^2));
end
Werr=Werr-diag(diag(Werr));
ISRef=Werr/N;
%SIRefica=-10*log10(sum(Werr,2)'/N);

% Werr=Werr/N+diag((ekurt-1)/(4*N));
% GCov=diag(Werr(:));
% for i=1:dim-1
%     for j=i+1:dim
%         Gjiij=-B(i)*B(j)/(J(i)^2*(B(j)+J(j)^2)+J(j)^2*B(i));
%         GCov((i-1)*dim+j,(j-1)*dim+i)=Gjiij;
%         GCov((j-1)*dim+i,(i-1)*dim+j)=Gjiij;
%     end
% end
Wefica=Wefica'*CC;

%varW_efica=reshape(diag(kron(Wefica',eye(dim))*GCov*kron(Wefica,eye(dim))),dim,dim);

%varA_efica=reshape(diag(kron(eye(dim),inv(Wefica))*GCov*kron(eye(dim),inv(Wefica'))),dim,dim);

function [g,gprime]=nonln(X)
[d N]=size(X);
g=zeros(d,N);
gprime=zeros(d,N);
h=1.06*N^(-1/5);
for i=1:d
    for k=1:N
        Z=(X(i,k)-X(i,:))/h;
        Zsquare=Z.^2;
        Phi=exp(-Zsquare/2)/sqrt(2*pi);
        f=mean(Phi,2)'/h;
        fprime=-mean(Z.*Phi,2)'/h^2;
        fprime2=mean((Zsquare-1).*Phi,2)'/h^3;
        g(i,k)=-fprime./f;
        gprime(i,k)=-fprime2./f+g(i,k)^2;
    end
end

function x=pwr(a,n)
x=a;
for i=2:n
    x=x.*a;
end