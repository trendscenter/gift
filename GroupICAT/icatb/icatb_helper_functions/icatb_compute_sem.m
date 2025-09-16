function sem = icatb_compute_sem(tmp_mean, tmp_ssq, N)
% Compute sem
sem = eps + ((sqrt((tmp_ssq - N*(tmp_mean.^2))/min([1, N]))) / sqrt(max([1, N - 1])));
sem = sem / sqrt(N);