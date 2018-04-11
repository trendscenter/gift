function [out, nm] = icatb_normit2(x1)
%normalizes x1 by L2 norm
%x1 = input vector
nm=norm(icatb_flatrow(x1(x1~=0)));
out = x1/nm;