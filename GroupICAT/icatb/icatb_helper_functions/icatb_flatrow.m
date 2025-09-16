function out = icatb_flatrow(in)
% out = icatb_flatrow(in)
%flatten into a row vector

out = reshape(in,prod(size(in)),1);
