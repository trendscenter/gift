function out = icatb_normdet(x1,p)
%detrends and normalizes x1
%x1 = input vector

sz = size(x1);
x1 = icatb_flatrow(x1);
if (~exist('p','var')),
   p=0;
end;

out = icatb_normit2(icatb_detrendit(x1,p));
out = reshape(out,sz);