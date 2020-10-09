function out = icatb_flatcol(in)
% out = icatb_flatcol
%flatten into a row vector

out = reshape(in,1,prod(size(in)));
