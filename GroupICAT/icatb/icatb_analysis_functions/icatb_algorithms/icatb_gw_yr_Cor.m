function gw=icatb_gw_yr_Cor(y,r,numsamp,threshold)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y is the output signal
% r is the reference signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
gw =-y*r'/numsamp-threshold;
