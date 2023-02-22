function pf = icatb_compute_var_evd(data, A, S)
% Function to compute percent variances using EVD method
%
% Inputs
% data - Voxels x subjects
% A - Subjects x components
% S - Voxels x components
%
% Outputs:
% pf - Percent variance
%

a = sum(eig(icatb_cov(data)));
vars = zeros(1, size(A, 2));
pf = vars;
for i = 1:size(A, 2)
    data_est = S(:, i)*A(:, i)';
    cov_m = icatb_cov(data_est);
    vars(i) = sum(eig(cov_m));
    res = data - data_est;
    cov_m = icatb_cov(res);
    var_res = sum(eig(cov_m)); %mean(sum(res.*res)/(size(res,2)-1));
    pf(i) = 100*(1 - (var_res/a));
end