function out = icatb_meanN(in)
% out = icatb_meanN(in)
% mean of image or matrix
N = length(size(in));

if (N == 2),
  out = mean(mean(in));
elseif (N == 3),
  out = mean(mean(mean(in)));
elseif (N == 4),
  out = mean(mean(mean(mean(in))));
end;
