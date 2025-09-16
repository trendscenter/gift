function M = icatb_compute_tvdfnc(FNC)
%% Compute temporal variation dFNC and stack dfnc and tvdfnc matrices
%
% Inputs:
% FNC - Windows x connectivity pairs
%
% Outputs:
% M - Stacked FNC and temporal variation dfnc
%

num_windows = size(FNC, 1);

tvdfnc = zeros(size(FNC));

% Forward difference at first timepoint
tvdfnc(1, :) = FNC(2, :) - FNC(1, :);

% central difference from window 2 to N - 1
for n = 2:num_windows - 1
    tvdfnc(n, :) = (FNC(n + 1, :) - FNC(n - 1, :)) / 2;
end

% backward difference at last window
tvdfnc(end, :) = FNC(end, :) - FNC(end-1, :);

FNC = diag(1./std(FNC,0,2))*FNC;
tvdfnc = diag(1./std(tvdfnc,0,2))*tvdfnc;

M = [FNC, tvdfnc];



