function c = initBicList()
%%%%
% To initialize the BicList global Structure 
%%%
measBuff = struct('subs',[],'comps',[],'freq',0); %initialize structure
global BicList;
BicList = measBuff;
end