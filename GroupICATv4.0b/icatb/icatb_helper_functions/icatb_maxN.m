function out = icatb_maxN(in)
%out = icatb_maxN(in)
%maximum of image or matrix
N = length(size(in));

if (N == 2),
  out = max(max(in));
elseif (N == 3),
  out = max(max(max(in)));
elseif (N == 4),
  out = max(max(max(max(in))));
end;
