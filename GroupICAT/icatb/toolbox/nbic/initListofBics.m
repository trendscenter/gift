function z = initListofBics()
%%%%
% To initialize the ListofBics global Structure 
%%%
measBuff = struct('subs',[],'comps',[],'freq',0); %initialize structure
global ListofBics;
ListofBics = measBuff;
end