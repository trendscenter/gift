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

function [ICOutMax, TCMax] = icatb_gigicar(FmriMatr, ICRefMax)
%written by Yuhui Du, CAS. 2012.

%Input:
%FmriMatr is the observed data with size of timepoints*voxelvolums£»
%ICRefMax includes the reference signals;

%Output£º
%ICOutMax includes the estimated ICs;
%TCMax is the obtained mixing matrix;


%thr = eps(class(FmriMatr)); %you can change the parameter. such as thr=0.02;

thr = eps(class(FmriMatr));
max_thr = 0.02;
stopping_tolerance = 1e-2;
thresholds = [thr, max_thr];

for ii = 1:length(thresholds)
    
    [ICOutMax, TCMax] = gig_icar(FmriMatr, ICRefMax, thresholds(ii));
    % make sure standard deviation of ic's to be 1
    std_ic = std(ICOutMax, 0, 2);
    err = norm(1 - std_ic);
    if (err > stopping_tolerance)
        disp('Standard deviation of the ic components are not unitary');
        if (ii == 1)
            disp(['Changing the min threshold of eigen values to ', num2str(thresholds(end)), ' and re-running gig-ica']);
        end
    else
        break;
    end
    
end

function [ICOutMax, TCMax, dpart] = gig_icar(FmriMatr, ICRefMax, thr)

%thr = eps(class(FmriMatr));

EsICnum = size(ICRefMax, 1);
[n, m] = size(FmriMatr);


meanDat = mean(FmriMatr,2);

FmriMatr = bsxfun(@minus, FmriMatr, meanDat);
%FmriMat = FmriMatr - repmat(mean(FmriMatr,2),[1,m]);
CovFmri = (FmriMatr*FmriMatr') / m;
[Esort, dsort] = eig(CovFmri);
dsort = diag(dsort);
filter_inds = find(sign(dsort) == -1);
dsort(filter_inds) = [];
Esort(:, filter_inds) = [];
dsort = abs(dsort);
%dsort = diag(dsort);
[dsort, flipped_inds] = sort(dsort, 'descend');
numpc = sum(dsort > thr);
Esort = Esort(:, flipped_inds);

Epart=Esort(:,1:numpc);
dpart=dsort(1:numpc);
Lambda_part=diag(dpart);
WhitenMatrix = (sqrtm(Lambda_part)) \ Epart';
Y=WhitenMatrix*FmriMatr;

if (thr < 1e-10 && numpc < n)
    Y = bsxfun(@rdivide, Y, std(Y, 0, 2));
end

ICRefMaxC = bsxfun(@minus, ICRefMax, mean(ICRefMax, 2));
ICRefMaxN = bsxfun(@rdivide, ICRefMaxC, std(ICRefMaxC, 0, 2));

NegeEva=zeros(EsICnum,1,class(FmriMatr));
for i=1:EsICnum
    NegeEva(i)=nege(ICRefMaxN(i,:));
end


YR = (1/m)*Y*ICRefMaxN';
referenceS = (ICRefMaxN*pinv(Y))';
referenceS = bsxfun(@rdivide, referenceS, sqrt(sum(referenceS.^2)));
iternum=100;
a=0.5;
b=1-a;
EGv=0.3745672075;
ErChuPai=2/pi;

WX = zeros(size(Y, 1), EsICnum, class(Y));


for ICnum=1:EsICnum
    reference=ICRefMaxN(ICnum,:);
    %wc=(reference*Yinv)';
    wc = referenceS(:, ICnum);
    %wc=wc/norm(wc);
    y1 = wc'*Y;
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
        tanhy1 = tanh(y1);
        EYgy=(1/m)*Y*(tanhy1)';
        %EYgy= sum(bsxfun(@times, Y,tanh(y1)/m),2);
        Jy1=(EGy1-EGv)^2;
        KwDaoshu=ErChuPai*c*(1/(1+(c*Jy1)^2));
        %Simgrad=(1/m)*Y*reference';
        Simgrad = YR(:,ICnum);
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
    WX(:, ICnum) = wx;
    %Source=wx'*Y;
    %ICOutMax(ICnum,:)=Source;
end

ICOutMax = WX'*Y;

FmriMatr = bsxfun(@plus, FmriMatr, meanDat);
TCMax = (1/m)*FmriMatr*ICOutMax';




function negentropy=nege(x)

y=log(cosh(x));
E1=mean(y);
E2=0.3745672075;
negentropy=(E1- E2)^2;
