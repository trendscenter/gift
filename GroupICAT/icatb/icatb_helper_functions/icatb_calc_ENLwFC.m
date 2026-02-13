function [ard_mat_subject] = icatb_calc_ENLwFC(data, s_option)
% This function can be used to calculate an explicitly nonlinear functional
% connectivity matrix for fMRI time series data. data: time*feature
% Cf. Kinsey et al. Networks extracted from nonlinear fMRI connectivity exhibit unique spatial variation and enhanced sensitivity to differences between individuals with schizophrenia and controls. Nat. Mental Health 2, 1464–1475 (2024). DOI: https://doi.org/10.1038/s44220-024-00341-y
% skinsey8@gsu.edu
% airaji@gsu.edu
% This code file is licensed under the MIT License (see LICENSE).

    
    
    switch lower(s_option)
        case 'nonlinear'
            b_linear = false;
        case 'linear'
            b_linear = true;
        otherwise
            error('icatb_calc_ENLwFC: only linear and non-linear option is supported.');
    end
    
    
    [nT,nF] = size(data);
    data = zscore(data);
    
    LINwFC = single(corr(data));
    LINwFC(isnan(LINwFC))=0;
    NLwFC = calc_dcorr(data);
    
    NLwFC = reshape(NLwFC,[nF*nF 1]);
    LINwFC = reshape(LINwFC,[nF*nF 1]);
    NLwFC = NLwFC - mean(NLwFC);
    LINwFC = LINwFC - mean(LINwFC);
    if b_linear
        ard_mat_subject = reshape(LINwFC,[nF nF]);
    else
        Beta = (LINwFC'*NLwFC)/(LINwFC'*LINwFC);
        ENLwFC = NLwFC - Beta*LINwFC; 
        ard_mat_subject = reshape(ENLwFC,[nF nF]);
    end
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