% *****************************************************************
% Copyright (c) Erik G. Learned-Miller, 2003.
% *****************************************************************
function [thetaStar,rotStar]=radicalOptTheta(x,stdev,m,reps,K,range)

% m is the number of intervals in an m-spacing
% reps is the number of points used in smoothing
% K is the number of angles theta to check for each Jacobi rotation.
[d,N]=size(x);

% This routine assumes that it gets whitened data.
% First, we augment the points with reps near copies of each point.

if reps==1
  xAug=x;
else
  xAug=randn(d,N*reps)*stdev+repmat(x,[1,reps]);
end

% Then we rotate this data to various angles, evaluate the sum of 
% the marginals, and take the min.

perc=range/(pi/2);
numberK=perc*K;
start=floor(K/2-numberK/2)+1;
endPt=ceil(K/2+numberK/2);

for i=1:K
  % Map theta from -pi/4 to pi/4 instead of 0 to pi/2.
  % This will allow us to use Amari-distance for test of
  % convergence.
  theta= (i-1)/(K-1)*pi/2-pi/4;
  rot=[cos(theta) -sin(theta); sin(theta) cos(theta)];
  rotPts=rot*xAug;

  for j=1:d
    marginalAtTheta(j) = icatb_vasicekm(rotPts(j,:),m);
  end
  ent(i)=sum(marginalAtTheta);
end

[val,ind]=sort(ent);
thetaStar= (ind(1)-1)/(K-1)*pi/2-pi/4;
fprintf(1,'rotated %5.2f degrees.\n',thetaStar/(2*pi)*360);
rotStar=[cos(thetaStar) -sin(thetaStar); sin(thetaStar) cos(thetaStar)];




