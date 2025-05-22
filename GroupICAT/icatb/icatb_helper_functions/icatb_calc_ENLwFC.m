function [ENLwFC] = icatb_calc_ENLwFC(data)
% This function can be used to calculate an explicitly nonlinear functional
% connectivity matrix for fMRI time series data. data: time*feature
% Cf. Kinsey et al. "Networks extracted from nonlinear fMRI connectivity
% exhibit unique spatial variation and enhanced sensitivity to differences
% between individuals with schizophrenia and controls" 
% skinsey8@gsu.edu
% airaji@gsu.edu
% This code file is licensed under the MIT License (see LICENSE).

[nT,nF] = size(data);
data = zscore(data);

LINwFC = single(corr(data));
NLwFC = calc_dcorr(data);

NLwFC = reshape(NLwFC,[nF*nF 1]);
LINwFC = reshape(LINwFC,[nF*nF 1]);
NLwFC = NLwFC - mean(NLwFC);
LINwFC = LINwFC - mean(LINwFC);
Beta = (LINwFC'*NLwFC)/(LINwFC'*LINwFC);
ENLwFC = NLwFC - Beta*LINwFC; 
ENLwFC = reshape(ENLwFC,[nF nF]);

end

function [Y] = calc_dcorr(X)
% This function can be used to calculate a distance correlation matrix for 
% fMRI time series data. X: time*feature

[nT,nF] = size(X);
Y = zeros(nT*nT,nF);

for ii = 1:nF
    x  = X(:,ii);
    a = pdist2(x,x);
    mcol = mean(a);
    A = a - mcol - mcol' + mean(mcol);
    Y(:,ii) = A(:);
end
clear X A

Y = single(Y);
Y = Y'*Y;
u = sqrt(diag(Y));
u = u*u';

Y = (Y./u);
Y(isnan(Y)) = 0;
clear u;
Y = sqrt(Y);

end