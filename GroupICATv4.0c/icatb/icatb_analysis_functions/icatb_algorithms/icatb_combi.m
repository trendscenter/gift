function [W, Wefica, Wwasobi, ISRwa, ISRef, metoda]= icatb_combi(X,AR_order,ortho)
%
% combines EFICA, WASOBI
%

if nargin<3
    ortho=false;
end

%COMMON PREPROCESSING
X = X-mean(X,2)*ones(1,size(X,2));
C = cov(X');
CC = C^(-1/2);
x = CC*X;


[Wefica, ISRef]=icatb_efica(x);
[Wwasobi, Wsobi,ISRwa]= icatb_ewasobi(x,10,0.9);
ISR1=sum(ISRef,2);
ISR2=sum(ISRwa,2);
imin=min([ISR1 ISR2]);
if imin(1)<imin(2)
    hotove=find(ISR1<imin(2));
    ostatni=find(ISR1>=imin(2));
    W=Wefica(hotove,:);
    if ortho
        Wost=null(Wefica(hotove,:))';
    else
        Wost=Wefica(ostatni,:);
    end
    y=Wost*x;
    metoda=ones(size(hotove'));
else
    hotove=find(ISR2<imin(1));
    ostatni=find(ISR2>=imin(1));
    W=Wwasobi(hotove,:);
    if ortho
        Wost=null(Wwasobi(hotove,:))';
    else
        Wost=Wwasobi(ostatni,:);
    end
    y=Wost*x;
    metoda=2*ones(size(hotove'));
end

while length(ostatni)>1
    [Wefi, ISRef0] = icatb_efica(y);
    [Wwa,AOL_init,ISRwa0]= icatb_ewasobi(y,10,0.9);
    ISR1=sum(ISRef0,2);
    ISR2=sum(ISRwa0,2);
    imin=min([ISR1 ISR2]);
    if imin(1)<imin(2)
        hotove=find(ISR1<imin(2));
        ostatni=find(ISR1>=imin(2));
        W=[W; Wefi(hotove,:)*Wost];
        if ortho
            Wost=null(Wefi(hotove,:))'*Wost;
        else
            Wost=Wefi(ostatni,:)*Wost;
        end
        y=Wost*x;
        metoda=[metoda ones(size(hotove'))];
    else
        hotove=find(ISR2<imin(1));
        ostatni=find(ISR2>=imin(1));
        W=[W; Wwa(hotove,:)*Wost];
        if ortho
            Wost=null(Wwa(hotove,:))'*Wost;
        else
            Wost=Wwa(ostatni,:)*Wost;
        end
        y=Wost*x;
        metoda=[metoda 2*ones(size(hotove'))];
    end
end
if ~isempty(ostatni)
    W=[W; Wost];
    metoda=[metoda 0];
end

W=W*CC;
Wefica=Wefica*CC;
Wwasobi=Wwasobi*CC;