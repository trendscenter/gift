function [name,a,b] = icatb_optimal5fun(data)
% [name,a,b] = icatb_optimal5fun(data)
% Created By: Baoming Hong, IOL, CT, USA
% Data: 1-08-04
%**************************************************************************
%*** the parameters: data denotes the estimated 1-D signal. f denotes the
%    estimated denisity function value, x is the corresponding axis values
%*** parameter b: distribution sacle factor parameter; a: mean/translation
%    factor, here not considering mean factor for Laplace, Logistic
%    Gaussian distribution, except for Logweilbull distribution...
%***                  Hong, Baoming, July 9, 2003
%*************************************************************************
% calcualte mean, variance, skewness & kurtosis
        PI = 3.1415926;
        meanS1 = mean(data');
        varS1 = var(data');
        skewS1 = skewness(data,0);
        kurS1 = kurtosis(data,0);
        b1 = sqrt(6*varS1)/3.1415926;
        a1 = meanS1 - b1*0.577215;
        
        meanS2 = mean(-data'); % for the signed data
        varS2 = var(-data');
        skewS2 = skewness(-data,0);
        kurS2 = kurtosis(-data,0);
        b2 = sqrt(6*varS2)/PI;
        a2 = meanS2 - b2*0.577215;
% now we want to estimate their distribution using kernel density estimation
       [f1, x1] = ksdensity(data);
       [f2, x2] = ksdensity(-data);   % note: the data: U=WX, probably be "signed"...
% initialization
        delt1 = zeros(4,1);
        delt2 = zeros(4,1);
        
% decide which distribution is optimal  
    for i=1:length(x1)  
         y1 = icatb_logistic_pdf(x1(i),0,1);  %%% here mean=0; var=1;
         y2 = icatb_lap_pdf(x1(i),meanS1,sqrt(varS1/2));
         y3 = normpdf(x1(i),meanS1,sqrt(varS1));
         y4 = icatb_logweibull_pdf(x1(i),a1,b1);
         delt1(1)=(y1-f1(i))^2 + delt1(1);
         delt1(2)=(y2-f1(i))^2 + delt1(2);
         delt1(3)=(y3-f1(i))^2 + delt1(3);
         delt1(4)=(y4-f1(i))^2 + delt1(4);
     end
         min1 = min(delt1);
     for i=1:length(x2)  
         y1 = icatb_logistic_pdf(x2(i),0,1);  %%% here mean=0; var=1;
         y2 = icatb_lap_pdf(x2(i),meanS2,sqrt(varS2/2));
         y3 = normpdf(x2(i),meanS2,sqrt(varS2));
         y4 = icatb_logweibull_pdf(x2(i),a2,b2);
         delt2(1)=(y1-f2(i))^2 + delt2(1);
         delt2(2)=(y2-f2(i))^2 + delt2(2);
         delt2(3)=(y3-f2(i))^2 + delt2(3);
         delt2(4)=(y4-f2(i))^2 + delt2(4);
     end
         min2 = min(delt2);
         
% Decide the minmum value and the corresponding position
    if    min1 < min2 | abs(min1-min2) <1.0e-6
          [mse,index] = min(delt1);
      else
          [mse,index] = min(delt2);
    end
      
% pick up the optimal function. If you expand or change the function, it is
% easy to update the codes below
    if index==1
        name = 'L';  %'logistic';
        a = meanS1;
        b = sqrt(3*varS1)/PI;
    elseif index==2
        name = 'P';  % 'laplace';
        a = meanS1;
        b = sqrt(varS1/2);
    elseif index ==3
        name = 'G';  %'gaussian';
        a = 0;
        b = 1;   % assigned the specific variance =1 here.
    elseif index==4
        name = 'W';  % 'logweibull';
        b = sqrt(6*varS1)/PI;
        a = meanS1 - b*0.577215;
    end