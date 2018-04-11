function [yfit] = icatb_getSplineFit(estimates,len,TR)

%TR = 2;
numP = floor(len/30);
t = 0:TR:(len-1)*TR;
t = t(:);
x0  = estimates;
yfit = x0(1)*t + x0(2)*t.^2;
for ii = 1:numP
    yfit = yfit + x0(2+ii) * sin(2*pi*ii*t/(len*TR)) + x0(2+ii) *cos(2*pi*ii*t/(len*TR));
end

%yfit = yfit';