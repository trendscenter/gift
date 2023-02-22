% *****************************************************************
% Copyright (c) Erik G. Learned-Miller, 2004.
% *****************************************************************
%
% Changes from version with release 1.0.
%
% 1) Switched to log from log2. log is slightly faster in Matlab.
% 2) Switched intvals calculation to make it more direct. This is a
%    substantial speed increase.

function h = icatb_vasicekm(v,m)

len=length(v);
vals=sort(v);

% Note that the intervals overlap for this estimator.
intvals=vals(m+1:len)-vals(1:len-m);
hvec=log(intvals);
h=sum(hvec);

