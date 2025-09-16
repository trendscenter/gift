function Gd=icatb_Gdfun(y) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--determining the first derivation of G(y)
%G(y)=exp(-y*y/2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %G(y)=exp(-a*y*y/2)
% a=0.001;
% Gd = -a*y.*exp(-a.*y.*y/2);

Gd = -y.*exp(-y.*y/2);