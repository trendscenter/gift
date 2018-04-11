% Contact ydu@mrn.org (yhdu@nlpr.ia.ac.cn) or yong.fan@ieee.org for bugs or questions 
%
%======================================================================================
%
%  Copyright (c) 2012 Yuhui DU and Yong FAN
%  All rights reserved.
%
% Redistribution and use in source or any other forms, with or without
% modification, are permitted provided that the following conditions are met:
%
%    * Redistributions of source code must retain the above copyright notice,
%      this list of conditions and the following disclaimer.
%
%    * Redistributions in any other form must reproduce the above copyright notice,
%      this list of conditions and the following disclaimer in the documentation
%      and/or other materials provided with the distribution.
%
%    * Neither the names of the copyright holders nor the names of future
%      contributors may be used to endorse or promote products derived from this
%      software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
% WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
% ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
% (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
% LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
% ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
% (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. 
%=====================================================================================

function [ICOutMax,TCMax] = icatb_gigicar(FmriMatr,ICRefMax)
%written by Yuhui Du, CAS. 2012.

%Input:
%FmriMatr is the observed data with size of timepoints*voxelvolums£»
%ICRefMax includes the reference signals;

%Output£º
%ICOutMax includes the estimated ICs;
%TCMax is the obtained mixing matrix;


[n,m]= size(FmriMatr);
FmriMat=FmriMatr - repmat(mean(FmriMatr,2),[1,m]);
CovFmri = (FmriMat*FmriMat') / m;
%CovFmri=cov(FmriMat',1);
[E,D]=eig(CovFmri);
EsICnum=size(ICRefMax,1); %EsICnum can be a number that is less than size(ICRefMax,1)
[eigenvalues index] = sort(diag(D));
cols=size(E,2);
Esort=zeros(size(E));
dsort=zeros(size(eigenvalues));
for i=1:cols
    Esort(:,i) = E(:, index(cols-i+1) );
    dsort(i)   = eigenvalues(index(cols-i+1) );
end

thr=0; %you can change the parameter. such as thr=0.02;
numpc=0;
for i=1:cols
    if dsort(i)>thr
        numpc=numpc+1;
    end
end

%dsum = sum(dsort);
% dsum_extract2 = sum(dsort(1:numpc));
% retained2=(dsum_extract2/dsum)*100;
% fprintf('%g%% of (non-zero) eigenvalues retained.\n', retained2);

Epart=Esort(:,1:numpc);
dpart=dsort(1:numpc);
Lambda_part=diag(dpart);
WhitenMatrix=(inv(sqrtm(Lambda_part)))*Epart';
Y=WhitenMatrix*FmriMat;

if thr<1e-10&&numpc<n
    for i=1:size(Y,1)
        Y(i,:)=Y(i,:)/std(Y(i,:));
    end
end

Yinv=pinv(Y);
ICRefMaxN=zeros(EsICnum,m);
ICRefMaxC=ICRefMax - repmat(mean(ICRefMax,2),[1,m]);
for i=1:EsICnum
    ICRefMaxN(i,:)=ICRefMaxC(i,:)/std(ICRefMaxC(i,:));
end

NegeEva=zeros(EsICnum,1);
for i=1:EsICnum
    NegeEva(i)=nege(ICRefMaxN(i,:));
end

iternum=100;
a=0.5;
b=1-a;
EGv=0.3745672075;
ErChuPai=2/pi;
ICOutMax=zeros(EsICnum,m);
for ICnum=1:EsICnum
    reference=ICRefMaxN(ICnum,:);
    wc=(reference*Yinv)';
    wc=wc/norm(wc);
    y1=wc'*Y;
    EyrInitial=(1/m)*(y1)*reference';
    NegeInitial=nege(y1);
    c=(tan((EyrInitial*pi)/2))/NegeInitial;
    IniObjValue=a*ErChuPai*atan(c*NegeInitial)+b*EyrInitial;
    
    itertime=1;
    Nemda=1;
    for i=1:iternum
        Cosy1=cosh(y1);
        logCosy1=log(Cosy1);
        EGy1=mean(logCosy1);
        Negama=EGy1-EGv;
        EYgy=(1/m)*Y*(tanh(y1))';
        Jy1=(EGy1-EGv)^2;
        KwDaoshu=ErChuPai*c*(1/(1+(c*Jy1)^2));
        Simgrad=(1/m)*Y*reference';
        g=a*KwDaoshu*2*Negama*EYgy+b*Simgrad;
        d=g/(g'*g)^0.5;
        wx=wc+Nemda*d;
        wx=wx/norm(wx);
        y3=wx'*Y;
        PreObjValue=a*ErChuPai*atan(c*nege(y3))+b*(1/m)*y3*reference';
        ObjValueChange=PreObjValue-IniObjValue;
        ftol=0.02;
        dg=g'*d;
        ArmiCondiThr=Nemda*ftol*dg;
        if ObjValueChange<ArmiCondiThr
            Nemda=Nemda/2;
            continue;
        end
        if (wc-wx)'*(wc-wx) <1.e-5
            break;
        else if itertime==iternum;
                break;
            end
        end
        IniObjValue=PreObjValue;
        y1=y3;
        wc=wx;
        itertime=itertime+1;
    end
    Source=wx'*Y;
    ICOutMax(ICnum,:)=Source;
end
TCMax=(1/m)*FmriMatr*ICOutMax';


function negentropy=nege(x)

y=log(cosh(x));
E1=mean(y);
E2=0.3745672075;
negentropy=(E1- E2)^2;