function out = icatb_detrendit(x1,p)
%function out = detrendit(x1)
%detrends x1 (removes mean)

  
if (~exist('p','var')),
   p=0;
end;

out = detrend(x1,p);