function [Z, Y] = icatb_ggmix(P, XDATA) 

XDATA = XDATA-P(8);
xpos = XDATA; xpos(find(XDATA < 0)) = 0;
xneg = XDATA; xneg(find(XDATA >= 0)) = 0;

%% gaussian

Y(1,:) = P(1)*exp(-XDATA.^2/(2*P(2)^2))/sqrt(2*pi*P(2)^2);


%% gamma 1
amp1 = gamma(P(4))*P(5)^P(4);
exponent1 = P(4)-1;

Y(2,:) = P(3)*(xpos.^exponent1).*exp(-xpos/P(5))/amp1;

%% gamma 2
amp2 = gamma(P(6))*P(7)^P(6);
exponent2 = P(6)-1;
Probgamma2 = max(1-P(1)-P(3), 0);
Y(3,:) = Probgamma2*((-xneg).^exponent2).*exp(xneg/P(7))/amp2;

%Y(3,:) = P(6)*((-xneg).^P(7)).*exp(xneg/P(8));

Z = sum(Y,1);