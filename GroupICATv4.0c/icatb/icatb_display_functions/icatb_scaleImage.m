function scaledData=icatb_scaleImage(tmin,tmax,imageData)
% scaledData=scaleImage(tmin,tmax,imageData)
% -------------------------------------------------------------------
% imagesData = (x*y*z,1)
% tmin = minValue to scale data to
% tmax = maxValue to scale data to


%scale data
% minVal = min(abs(imageData(find(imageData~=0))));
% maxVal = max(abs(imageData(find(imageData~=0))));
minVal = min(imageData);
maxVal = max(imageData);
rangeVal = maxVal-minVal;
trange = tmax-tmin;
if (rangeVal == 0)
    rangeVal = eps;
end
scaledData = imageData;
scaledData = (((imageData-minVal)./rangeVal)./(1/trange))+tmin;
%scaledData(find(imageData~=0)) = (((abs(imageData(find(imageData~=0)))-minVal)./rangeVal)./(1/trange))+tmin;


