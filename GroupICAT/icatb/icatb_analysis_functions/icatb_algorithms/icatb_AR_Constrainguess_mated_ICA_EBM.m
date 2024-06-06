function [W, rho, Cost] = icatb_AR_Constrainguess_mated_ICA_EBM(X, guess_mat, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adaptive Reverse Constrained EBM
% ICA-EBM: ICA by Entropy Bound Minimization (real-valued version)
% Four nonlinearities
% x^4,  |x|/(1+|x|),    x|x|/(10+|x|),  and     x/(1+x^2)
% are used for entropy bound calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs:
% X:    mixtures
% guess_mat: matrix where each column is a reference vector. num_guess
%       is the number of reference signals.
% output:
% W:    demixing matrix
%
% rho: tuned constrained parameter at each iteration (num_guess x K x iter)
% matrix where num_guess is number of constraints, iter is number of
% iteration algorithm runs
%
% Cost: number of iterations it took to get results
%
% Program by Hanlu Yang. Please contact me at hyang3@umbc.edu
% Example:  [W, rho_n_arr, ~]= icatb_AR_Constrained_ICA_EBM(X,guess_mat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Modification (6/4/2024 Cyrus) to flip signs depending on reference 
% components for improved correction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:

% [1] Yang, H., Ghayem F., Gabrielson B., M. A. B. S. Akhonda, Calhoun, V. D., and  Adali, T.,
% "Constrained independent component analysis based on entropy bound minimization 
% for subgroup identification from multisubject fMRI data"

% [2] Vu, T., Laport, F., Yang, H., and Adali, T., 
% "Constrained Independent Vector Analysis with References: 
% Algorithms and Performance Evaluation." In 2023 57th Asilomar Conference on Signals, 
% Systems, and Computers (pp. 827-831). IEEE.

% [3] Yang, H., Ghayem F., Gabrielson B., M. A. B. S. Akhonda, Calhoun, V. D., and  Adali, T.,
% "Constrained independent component analysis based on entropy bound minimization
% for subgroup identification from multisubject fMRI data." 
% In IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP). IEEE, 2023.

% [4] Li, X. L., & Adali, T. (2009, September). A novel entropy estimator and its application to ICA.
% In 2009 IEEE International Workshop on Machine Learning for Signal Processing (pp. 1-6). IEEE.
%
% [5] Li, X. L., & Adali, T. (2010). Independent component analysis by entropy bound minimization.
% IEEE Transactions on Signal Processing, 58(10), 5151-5164.
%
%
%

xOrig=X;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Part 0: Here begins some pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build default options structure
options = struct( ...
    'opt_approach', 'newton', ... % optimization type: gradient, (newton), quasi
    'complex_valued', false, ... % if any input is complex or quasi approach used then setting is forced to true
    'circular', false, ... % set to true to only consider circular for complex-valued cost function
    'whiten', true, ... % whitening is optional (except for quasi approach it is required)
    'verbose', false, ... % verbose true enables print statements
    'A', [], ... % true mixing matrices A, automatically sets verbose
    'W_init', [], ... % initial estimates for demixing matrices in W
    'jdiag_initW', false, ... % use CCA (K=2) or joint diagonalization (K>2)
    'maxIter', 512, ... % max number of iterations
    'WDiffStop', 1e-6, ... % stopping criterion
    'alpha0', 1.0, ... % initial step size scaling (will be doubled for complex-valued)
    'Save_W', false, ...
    'gam', 100 ... % lagrange multiplier step parameter (default 3)F
    );

% load in user supplied options
options = getopt(options, varargin{:});

supplyA = ~isempty(options.A); % set to true if user has supplied true mixing matrices
blowup = 1e3;
alphaScale = 0.9; % alpha0 to alpha0*alphaScale when cost does not decrease
%options.rho = rho;
alphaMin = options.WDiffStop; % alpha0 will max(alphaMin,alphaScale*alpha0)
outputISI = false;
opt_approach = find(strcmp(options.opt_approach, ...
    {'gradient', 'newton', 'quasi'}), 1); % 1='gradient', 2='newton', 3='quasi'
% options.whiten=(options.whiten || (opt_approach==3)); % whitening required for real-valued quasi & complex-valued gradient approach

rho = 0:0.01:1; % constraint parameter rho is ranging from 0.1 to 0.99 with 0.01 step


max_iter_fastica = 100;
max_iter_orth = 1000;
max_iter_orth_refine = 1000;
max_iter_nonorth = 1000;
saddle_test_enable = 1;
tolerance = 1e-6;
max_cost_increase_number = 5;
stochastic_search_factor = 1;
%mu_c = (1-median(rho))*options.gam ;  % initialize mu_c with median of rho

gam = 100;

verbose = 0; % report the progress if verbose==1
show_cost = 0; % show the cost values vs. iterations at each stage if show_cost==1

% Load 8 measuring functions. But we only use 4 of them.
K = 8;
load icatb_nf_table;

[N, T] = size(X);
[X, P] = pre_processing(X);

%% Initialize W
if ~isempty(options.W_init)
    W = options.W_init;
else

    % use SEA, or fastica with pow3, to provide an initial guess
    W = randn(N, N);
end

W = symdecor(W);
num_guess = size(guess_mat, 2);
mu_c = max(0, randn(num_guess, 1));
%mu_c = mu_c * ones(num_guess,1);
num_W = size(W, 1);
%Re-sort existing W based on correlation with matrix W:
for kl = 1:num_guess
    r_n_c = guess_mat(:, kl);
    for lp = 1:num_W
        w = W(lp, :).';
        R_par = corrcoef(X'*w, r_n_c); %corr(y',rnc)
        corr_w_guess(kl, lp) = R_par(1, 2);
    end
end
%We may need to use auction to chose order:
[~, max_idx] = max(abs(corr_w_guess), [], 2); % ; Added by Zois

if length(unique(max_idx)) ~= num_guess
    [colsol, ~] = auction((1 - abs(corr_w_guess))');
    max_idx = colsol';
end

c = setxor((1:num_W).', max_idx);
sort_order = [max_idx; c];

W = W(sort_order, :);

last_W = W;
best_W = W;
Cost = zeros(max_iter_fastica, 1);
min_cost = inf;
cost_increase_counter = 0;
for iter = 1:max_iter_fastica

    Y = W * X;

    for n = 1:N

        y = Y(n, :); %y = w'*X;

        % evaluate the upper bound of negentropy of the nth component
        NE_Bound = zeros(K, 1);
        EGx = zeros(K, 1);

        % we only need to calculate these quantities once
        yy = y .* y;
        sign_y = sign(y);
        abs_y = sign_y .* y;
        inv_pabs_y = 1 ./ (1 + abs_y);
        inv_pabs_yy = 1 ./ (1 + yy);
        inv_p10abs_y = 1 ./ (10 + abs_y);

        % G1(x) = x^4
        EGx(1) = sum(yy.*yy) / T;
        if EGx(1) < nf1.min_EGx
            NE_Bound(1) = simplified_ppval(nf1.pp_slope, nf1.min_EGx) * (EGx(1) - nf1.min_EGx);
            NE_Bound(1) = simplified_ppval(nf1.pp, nf1.min_EGx) + abs(NE_Bound(1));
        else if EGx(1) > nf1.max_EGx
                % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                NE_Bound(1) = 0;
        else
            NE_Bound(1) = simplified_ppval(nf1.pp, EGx(1));
        end
        end

        % G3 = @(x) abs(x)./( 1+abs(x) );
        EGx(3) = 1 - sum(inv_pabs_y) / T;
        if EGx(3) < nf3.min_EGx
            NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.min_EGx) * (EGx(3) - nf3.min_EGx);
            NE_Bound(3) = simplified_ppval(nf3.pp, nf3.min_EGx) + abs(NE_Bound(3));
        else if EGx(3) > nf3.max_EGx
                NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.max_EGx) * (EGx(3) - nf3.max_EGx);
                NE_Bound(3) = simplified_ppval(nf3.pp, nf3.max_EGx) + abs(NE_Bound(3));
        else
            NE_Bound(3) = simplified_ppval(nf3.pp, EGx(3));
        end
        end

        % G5 = @(x) x*|x|/(10+|x|);
        EGx(5) = sum(y.*abs_y.*inv_p10abs_y) / T;
        if EGx(5) < nf5.min_EGx
            NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.min_EGx) * (EGx(5) - nf5.min_EGx);
            NE_Bound(5) = simplified_ppval(nf5.pp, nf5.min_EGx) + abs(NE_Bound(5));
        else if EGx(5) > nf5.max_EGx
                NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.max_EGx) * (EGx(5) - nf5.max_EGx);
                NE_Bound(5) = simplified_ppval(nf5.pp, nf5.max_EGx) + abs(NE_Bound(5));
        else
            NE_Bound(5) = simplified_ppval(nf5.pp, EGx(5));
        end
        end

        % G7 = @(x) x/(1+x^2);
        EGx(7) = sum(y.*inv_pabs_yy) / T;
        if EGx(7) < nf7.min_EGx
            NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.min_EGx) * (EGx(7) - nf7.min_EGx);
            NE_Bound(7) = simplified_ppval(nf7.pp, nf7.min_EGx) + abs(NE_Bound(7));
        else if EGx(7) > nf7.max_EGx
                NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.max_EGx) * (EGx(7) - nf7.max_EGx);
                NE_Bound(7) = simplified_ppval(nf7.pp, nf7.max_EGx) + abs(NE_Bound(7));
        else
            NE_Bound(7) = simplified_ppval(nf7.pp, EGx(7));
        end
        end

        % select the tightest bound
        [max_NE, max_i] = max(NE_Bound);
        negentropy_array(n) = max_NE;
        Cost(iter) = Cost(iter) - max_NE;

    end


    if Cost(iter) < min_cost
        min_cost = Cost(iter);
        best_W = last_W;
        cost_increase_counter = 0;
    else
        cost_increase_counter = cost_increase_counter + 1;
    end

    W = (Y .* Y .* Y) * X' / T - 3 * W;
    W = symdecor(W);

    if 1 - min(abs(diag(W*last_W'))) < tolerance
        break;
    else
        last_W = W;
    end

    if cost_increase_counter > max_cost_increase_number
        break;
    end

end
W = best_W;
if show_cost
    figure;
    subplot(1, 4, 1);
    plot(1:iter, Cost(1:iter));
    xlabel('Number of iterations')
    ylabel('Cost')
    title('(0) Initialize use SEA')
end

% [End of Part 0]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Part I: orthogonal ICA
%  varying step size, stochastic gradient search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    fprintf('\nOrthogonal ICA stage.');
end

num_guess = size(guess_mat, 2);
num_W = size(W, 1);
%Re-sort existing W based on correlation with matrix W:
for kl = 1:num_guess
    r_n_c = guess_mat(:, kl);
    for lp = 1:num_W
        w = W(lp, :).';
        R_par = corrcoef(X'*w, r_n_c); %corr(y',rnc)
        corr_w_guess(kl, lp) = R_par(1, 2);
    end
end
%We may need to use auction to chose order:
[~, max_idx] = max(abs(corr_w_guess), [], 2); % ; Added by Zois

if length(unique(max_idx)) ~= num_guess
    [colsol, ~] = auction((1 - abs(corr_w_guess))');
    max_idx = colsol';
end

c = setxor((1:num_W).', max_idx);
sort_order = [max_idx; c];

W = W(sort_order, :);

last_W = W;
best_W = W;
Cost = zeros(max_iter_orth, 1);
min_cost = inf;
min_cost_queue = min_cost * ones(max_iter_orth, 1);
mu = 1 / 6.25;
min_mu = 1 / 50;
cost_increase_counter = 0;
fastica_on = 1;
error = 0;
max_negentropy = zeros(N, 1);
negentropy_array = zeros(N, 1);
for iter = 1:max_iter_orth

    Y = W * X;
    for n = 1:N

        w = W(n, :)';
        y = Y(n, :); %y = w'*X;

        % evaluate the upper bound of negentropy of the nth component
        NE_Bound = zeros(K, 1);
        EGx = zeros(K, 1);

        % we only need to calculate these quantities once
        yy = y .* y;
        sign_y = sign(y);
        abs_y = sign_y .* y;
        inv_pabs_y = 1 ./ (1 + abs_y);
        inv_pabs_yy = 1 ./ (1 + yy);
        inv_p10abs_y = 1 ./ (10 + abs_y);

        % G1(x) = x^4
        EGx(1) = sum(yy.*yy) / T;
        if EGx(1) < nf1.min_EGx
            NE_Bound(1) = simplified_ppval(nf1.pp_slope, nf1.min_EGx) * (EGx(1) - nf1.min_EGx);
            NE_Bound(1) = simplified_ppval(nf1.pp, nf1.min_EGx) + abs(NE_Bound(1));
        else if EGx(1) > nf1.max_EGx
                % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                NE_Bound(1) = 0;
        else
            NE_Bound(1) = simplified_ppval(nf1.pp, EGx(1));
        end
        end

        % G3 = @(x) abs(x)./( 1+abs(x) );
        EGx(3) = 1 - sum(inv_pabs_y) / T;
        if EGx(3) < nf3.min_EGx
            NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.min_EGx) * (EGx(3) - nf3.min_EGx);
            NE_Bound(3) = simplified_ppval(nf3.pp, nf3.min_EGx) + abs(NE_Bound(3));
        else if EGx(3) > nf3.max_EGx
                NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.max_EGx) * (EGx(3) - nf3.max_EGx);
                NE_Bound(3) = simplified_ppval(nf3.pp, nf3.max_EGx) + abs(NE_Bound(3));
        else
            NE_Bound(3) = simplified_ppval(nf3.pp, EGx(3));
        end
        end

        % G5 = @(x) x*|x|/(10+|x|);
        EGx(5) = sum(y.*abs_y.*inv_p10abs_y) / T;
        if EGx(5) < nf5.min_EGx
            NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.min_EGx) * (EGx(5) - nf5.min_EGx);
            NE_Bound(5) = simplified_ppval(nf5.pp, nf5.min_EGx) + abs(NE_Bound(5));
        else if EGx(5) > nf5.max_EGx
                NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.max_EGx) * (EGx(5) - nf5.max_EGx);
                NE_Bound(5) = simplified_ppval(nf5.pp, nf5.max_EGx) + abs(NE_Bound(5));
        else
            NE_Bound(5) = simplified_ppval(nf5.pp, EGx(5));
        end
        end

        % G7 = @(x) x/(1+x^2);
        EGx(7) = sum(y.*inv_pabs_yy) / T;
        if EGx(7) < nf7.min_EGx
            NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.min_EGx) * (EGx(7) - nf7.min_EGx);
            NE_Bound(7) = simplified_ppval(nf7.pp, nf7.min_EGx) + abs(NE_Bound(7));
        else if EGx(7) > nf7.max_EGx
                NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.max_EGx) * (EGx(7) - nf7.max_EGx);
                NE_Bound(7) = simplified_ppval(nf7.pp, nf7.max_EGx) + abs(NE_Bound(7));
        else
            NE_Bound(7) = simplified_ppval(nf7.pp, EGx(7));
        end
        end

        % select the tightest bound
        [max_NE, max_i] = max(NE_Bound);
        negentropy_array(n) = max_NE;
        Cost(iter) = Cost(iter) - max_NE;

        if ~fastica_on
            weight = rand(1, T);
        end

        % Perform orthogonal ICA
        switch max_i

            case 1

                % G1(x) = x^4
                % g1(x) = 4x^3
                % g1' =  12x^2
                if fastica_on
                    grad = X * (4 * y .* yy)' / T;
                    Edgx = 12;
                else
                    grad = X * (4 * weight .* y .* yy)' / sum(weight);
                    vEGx = 2 * (EGx(1) > nf1.critical_point) - 1;
                end

            case 3

                % G3 = @(x) abs(x)./( 1+abs(x) );
                % g3 = @(x) sign(x)./( 1+abs(x) ).^2;
                % g3' = -2/(1+abs(x))^3
                if fastica_on
                    grad = X * (sign_y .* inv_pabs_y .* inv_pabs_y)' / T; %grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
                    Edgx = sum(-2*inv_pabs_y.*inv_pabs_y.*inv_pabs_y) / T; % Edgx = sum(-2*inv_pabs_y.^3)/T;
                else
                    grad = X * (weight .* sign_y .* inv_pabs_y .* inv_pabs_y)' / sum(weight); %grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
                    vEGx = 2 * (EGx(3) > nf3.critical_point) - 1;
                end

            case 5

                % G5 = @(x) x*|x|/(10+|x|);
                % g5 = @(x) |x|*(20+|x|)/(10+|x|)^2;
                % g5 = 200*sign(x) / (10+|x|)^3
                if fastica_on
                    grad = X * (abs_y .* (20 + abs_y) .* inv_p10abs_y.^2)' / T;
                    Edgx = sum(200*sign_y.*inv_p10abs_y.*inv_p10abs_y.*inv_p10abs_y) / T; %Edgx = sum( 200*sign_y.*inv_p10abs_y.^3 )/T;
                else
                    grad = X * (weight .* abs_y .* (20 + abs_y) .* inv_p10abs_y.^2)' / sum(weight);
                    vEGx = 2 * (EGx(5) > nf5.critical_point) - 1;
                end

            case 7

                % G7 = @(x) x/(1+x^2)
                % g7 = @(x) (1-x^2)/(1+x^2)^2
                % g7' = 2x(x^2-3)/(1+x^2)^3
                if fastica_on
                    grad = X * ((1 - yy) .* inv_pabs_yy.^2)' / T;
                    Edgx = sum(2*y.*(yy - 3).*inv_pabs_yy.*inv_pabs_yy.*inv_pabs_yy) / T; %Edgx = sum( 2*y.*(yy-3).*inv_pabs_yy.^3 )/T;
                else
                    grad = X * (weight .* (1 - yy) .* inv_pabs_yy.^2)' / sum(weight);
                    vEGx = 2 * (EGx(7) > nf7.critical_point) - 1;
                end

            otherwise
                ;
        end


        if fastica_on
            w1 = grad - Edgx * w;
        else
            grad = vEGx * grad;
            grad = grad - (w' * grad) * w;
            grad = grad / norm(grad);
            w1 = w + mu * grad;
        end

        W(n, :) = w1';

    end

    W = symdecor(W);

    if Cost(iter) < min_cost
        cost_increase_counter = 0;
        min_cost = Cost(iter);
        best_W = last_W;
        max_negentropy = negentropy_array;
    else
        cost_increase_counter = cost_increase_counter + 1;
    end
    min_cost_queue(iter) = min_cost;

    if fastica_on
        if cost_increase_counter >= max_cost_increase_number || 1 - min(abs(diag(W*last_W'))) < tolerance
            cost_increase_counter = 0;
            W = best_W;
            last_W = W;
            iter_fastica = iter;
            fastica_on = 0;
            continue;
        end
    else
        if cost_increase_counter > stochastic_search_factor * max_cost_increase_number
            if mu > min_mu
                cost_increase_counter = 0;
                W = best_W;
                last_W = W;
                mu = mu / 2;
                continue;
            else
                break;
            end
        end
    end

    last_W = W;

end
% [End of Part I]
W = best_W;

if show_cost
    subplot(1, 4, 2);
    plot([1:iter_fastica], Cost([1:iter_fastica]), 'r');
    hold on;
    plot([iter_fastica + 1:iter], Cost([iter_fastica + 1:iter]), 'b');
    xlabel('number of iterations')
    ylabel('Cost')
    title('(a) Orthogonal ICA')
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           PartII: check saddle points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if saddle_test_enable
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose
        fprintf('\nSaddle point detection.');
    end
    SADDLE_TESTED = 0;
    saddle_tested = 1;
    while saddle_tested

        saddle_tested = 0;
        Y = W * X;
        for m = 1:N
            w1 = W(m, :)';
            ym = Y(m, :); %ym = w1'*X;
            for n = m + 1:N
                w2 = W(n, :)';
                yn = Y(n, :); %yn = w2'*X;

                yr1 = (ym + yn) / sqrt(2);
                yr2 = (ym - yn) / sqrt(2);

                y = yr1;
                % we only need to calculate these quantities once
                yy = y .* y;
                sign_y = sign(y);
                abs_y = sign_y .* y;
                inv_pabs_y = 1 ./ (1 + abs_y);
                inv_pabs_yy = 1 ./ (1 + yy);
                inv_p10abs_y = 1 ./ (10 + abs_y);

                % G1(x) = x^4
                EGx(1) = sum(yy.*yy) / T;
                if EGx(1) < nf1.min_EGx
                    NE_Bound(1) = simplified_ppval(nf1.pp_slope, nf1.min_EGx) * (EGx(1) - nf1.min_EGx);
                    NE_Bound(1) = simplified_ppval(nf1.pp, nf1.min_EGx) + abs(NE_Bound(1));
                else if EGx(1) > nf1.max_EGx
                        % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                        % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                        NE_Bound(1) = 0;
                else
                    NE_Bound(1) = simplified_ppval(nf1.pp, EGx(1));
                end
                end

                % G3 = @(x) abs(x)./( 1+abs(x) );
                EGx(3) = 1 - sum(inv_pabs_y) / T;
                if EGx(3) < nf3.min_EGx
                    NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.min_EGx) * (EGx(3) - nf3.min_EGx);
                    NE_Bound(3) = simplified_ppval(nf3.pp, nf3.min_EGx) + abs(NE_Bound(3));
                else if EGx(3) > nf3.max_EGx
                        NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.max_EGx) * (EGx(3) - nf3.max_EGx);
                        NE_Bound(3) = simplified_ppval(nf3.pp, nf3.max_EGx) + abs(NE_Bound(3));
                else
                    NE_Bound(3) = simplified_ppval(nf3.pp, EGx(3));
                end
                end

                % G5 = @(x) x*|x|/(10+|x|);
                EGx(5) = sum(y.*abs_y.*inv_p10abs_y) / T;
                if EGx(5) < nf5.min_EGx
                    NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.min_EGx) * (EGx(5) - nf5.min_EGx);
                    NE_Bound(5) = simplified_ppval(nf5.pp, nf5.min_EGx) + abs(NE_Bound(5));
                else if EGx(5) > nf5.max_EGx
                        NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.max_EGx) * (EGx(5) - nf5.max_EGx);
                        NE_Bound(5) = simplified_ppval(nf5.pp, nf5.max_EGx) + abs(NE_Bound(5));
                else
                    NE_Bound(5) = simplified_ppval(nf5.pp, EGx(5));
                end
                end

                % G7 = @(x) x/(1+x^2);
                EGx(7) = sum(y.*inv_pabs_yy) / T;
                if EGx(7) < nf7.min_EGx
                    NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.min_EGx) * (EGx(7) - nf7.min_EGx);
                    NE_Bound(7) = simplified_ppval(nf7.pp, nf7.min_EGx) + abs(NE_Bound(7));
                else if EGx(7) > nf7.max_EGx
                        NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.max_EGx) * (EGx(7) - nf7.max_EGx);
                        NE_Bound(7) = simplified_ppval(nf7.pp, nf7.max_EGx) + abs(NE_Bound(7));
                else
                    NE_Bound(7) = simplified_ppval(nf7.pp, EGx(7));
                end
                end

                % select the tightest bound
                [max_NE, max_i] = max(NE_Bound);
                negentropy1 = max_NE;

                y = yr2;
                % we only need to calculate these quantities once
                yy = y .* y;
                sign_y = sign(y);
                abs_y = sign_y .* y;
                inv_pabs_y = 1 ./ (1 + abs_y);
                inv_pabs_yy = 1 ./ (1 + yy);
                inv_p10abs_y = 1 ./ (10 + abs_y);

                % G1(x) = x^4
                EGx(1) = sum(yy.*yy) / T;
                if EGx(1) < nf1.min_EGx
                    NE_Bound(1) = simplified_ppval(nf1.pp_slope, nf1.min_EGx) * (EGx(1) - nf1.min_EGx);
                    NE_Bound(1) = simplified_ppval(nf1.pp, nf1.min_EGx) + abs(NE_Bound(1));
                else if EGx(1) > nf1.max_EGx
                        % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                        % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                        NE_Bound(1) = 0;
                else
                    NE_Bound(1) = simplified_ppval(nf1.pp, EGx(1));
                end
                end

                % G3 = @(x) abs(x)./( 1+abs(x) );
                EGx(3) = 1 - sum(inv_pabs_y) / T;
                if EGx(3) < nf3.min_EGx
                    NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.min_EGx) * (EGx(3) - nf3.min_EGx);
                    NE_Bound(3) = simplified_ppval(nf3.pp, nf3.min_EGx) + abs(NE_Bound(3));
                else if EGx(3) > nf3.max_EGx
                        NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.max_EGx) * (EGx(3) - nf3.max_EGx);
                        NE_Bound(3) = simplified_ppval(nf3.pp, nf3.max_EGx) + abs(NE_Bound(3));
                else
                    NE_Bound(3) = simplified_ppval(nf3.pp, EGx(3));
                end
                end

                % G5 = @(x) x*|x|/(10+|x|);
                EGx(5) = sum(y.*abs_y.*inv_p10abs_y) / T;
                if EGx(5) < nf5.min_EGx
                    NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.min_EGx) * (EGx(5) - nf5.min_EGx);
                    NE_Bound(5) = simplified_ppval(nf5.pp, nf5.min_EGx) + abs(NE_Bound(5));
                else if EGx(5) > nf5.max_EGx
                        NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.max_EGx) * (EGx(5) - nf5.max_EGx);
                        NE_Bound(5) = simplified_ppval(nf5.pp, nf5.max_EGx) + abs(NE_Bound(5));
                else
                    NE_Bound(5) = simplified_ppval(nf5.pp, EGx(5));
                end
                end

                % G7 = @(x) x/(1+x^2);
                EGx(7) = sum(y.*inv_pabs_yy) / T;
                if EGx(7) < nf7.min_EGx
                    NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.min_EGx) * (EGx(7) - nf7.min_EGx);
                    NE_Bound(7) = simplified_ppval(nf7.pp, nf7.min_EGx) + abs(NE_Bound(7));
                else if EGx(7) > nf7.max_EGx
                        NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.max_EGx) * (EGx(7) - nf7.max_EGx);
                        NE_Bound(7) = simplified_ppval(nf7.pp, nf7.max_EGx) + abs(NE_Bound(7));
                else
                    NE_Bound(7) = simplified_ppval(nf7.pp, EGx(7));
                end
                end

                % select the tightest bound
                [max_NE, max_i] = max(NE_Bound);
                negentropy2 = max_NE;

                if negentropy1 + negentropy2 > max_negentropy(m) + max_negentropy(n) + eps
                    if verbose
                        fprintf('\nRotating %g %g', m, n);
                    end
                    max_negentropy(m) = negentropy1;
                    max_negentropy(n) = negentropy2;
                    W(m, :) = (w1 + w2)' / sqrt(2);
                    W(n, :) = (w1 - w2)' / sqrt(2);
                    Y(m, :) = yr1;
                    Y(n, :) = yr2;
                    ym = yr1;
                    w1 = W(m, :)';
                    saddle_tested = 1;
                    SADDLE_TESTED = 1;
                end

            end
        end

    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    SADDLE_TESTED = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if SADDLE_TESTED == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %           Part III: if saddles are detected, refine orthogonal ICA
    %  fixed step size, gradient search
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verbose
        fprintf('\nOrthogonal ICA refinement is required because saddles are detected.');
    end
    last_W = W;
    best_W = W;
    Cost = zeros(max_iter_orth_refine, 1);
    min_cost = inf;
    min_cost_queue = min_cost * ones(max_iter_orth_refine, 1);
    mu = 1 / 50;
    cost_increase_counter = 0;
    fastica_on = 1;
    error = 0;
    for iter = 1:max_iter_orth_refine

        for n = 1:N

            w = W(n, :)';
            y = w' * X;

            % evaluate the upper bound of negentropy of the nth component
            NE_Bound = zeros(K, 1);
            EGx = zeros(K, 1);

            % we only need to calculate these quantities once
            yy = y .* y;
            sign_y = sign(y);
            abs_y = sign_y .* y;
            inv_pabs_y = 1 ./ (1 + abs_y);
            inv_pabs_yy = 1 ./ (1 + yy);
            inv_p10abs_y = 1 ./ (10 + abs_y);

            % G1(x) = x^4
            EGx(1) = sum(yy.*yy) / T;
            if EGx(1) < nf1.min_EGx
                NE_Bound(1) = simplified_ppval(nf1.pp_slope, nf1.min_EGx) * (EGx(1) - nf1.min_EGx);
                NE_Bound(1) = simplified_ppval(nf1.pp, nf1.min_EGx) + abs(NE_Bound(1));
            else if EGx(1) > nf1.max_EGx
                    % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                    % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                    NE_Bound(1) = 0;
            else
                NE_Bound(1) = simplified_ppval(nf1.pp, EGx(1));
            end
            end

            % G3 = @(x) abs(x)./( 1+abs(x) );
            EGx(3) = 1 - sum(inv_pabs_y) / T;
            if EGx(3) < nf3.min_EGx
                NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.min_EGx) * (EGx(3) - nf3.min_EGx);
                NE_Bound(3) = simplified_ppval(nf3.pp, nf3.min_EGx) + abs(NE_Bound(3));
            else if EGx(3) > nf3.max_EGx
                    NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.max_EGx) * (EGx(3) - nf3.max_EGx);
                    NE_Bound(3) = simplified_ppval(nf3.pp, nf3.max_EGx) + abs(NE_Bound(3));
            else
                NE_Bound(3) = simplified_ppval(nf3.pp, EGx(3));
            end
            end

            % G5 = @(x) x*|x|/(10+|x|);
            EGx(5) = sum(y.*abs_y.*inv_p10abs_y) / T;
            if EGx(5) < nf5.min_EGx
                NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.min_EGx) * (EGx(5) - nf5.min_EGx);
                NE_Bound(5) = simplified_ppval(nf5.pp, nf5.min_EGx) + abs(NE_Bound(5));
            else if EGx(5) > nf5.max_EGx
                    NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.max_EGx) * (EGx(5) - nf5.max_EGx);
                    NE_Bound(5) = simplified_ppval(nf5.pp, nf5.max_EGx) + abs(NE_Bound(5));
            else
                NE_Bound(5) = simplified_ppval(nf5.pp, EGx(5));
            end
            end

            % G7 = @(x) x/(1+x^2);
            EGx(7) = sum(y.*inv_pabs_yy) / T;
            if EGx(7) < nf7.min_EGx
                NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.min_EGx) * (EGx(7) - nf7.min_EGx);
                NE_Bound(7) = simplified_ppval(nf7.pp, nf7.min_EGx) + abs(NE_Bound(7));
            else if EGx(7) > nf7.max_EGx
                    NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.max_EGx) * (EGx(7) - nf7.max_EGx);
                    NE_Bound(7) = simplified_ppval(nf7.pp, nf7.max_EGx) + abs(NE_Bound(7));
            else
                NE_Bound(7) = simplified_ppval(nf7.pp, EGx(7));
            end
            end

            % select the tightest bound
            [max_NE, max_i] = max(NE_Bound);
            negentropy_array(n) = max_NE;
            Cost(iter) = Cost(iter) - max_NE;

            % Perform orthogonal ICA
            switch max_i

                case 1

                    % G1(x) = x^4
                    % g1(x) = 4x^3
                    % g1' =  12x^2
                    grad = X * (4 * y .* yy)' / T;
                    if fastica_on
                        Edgx = 12;
                    else
                        vEGx = 2 * (EGx(1) > nf1.critical_point) - 1;
                    end

                case 3

                    % G3 = @(x) abs(x)./( 1+abs(x) );
                    % g3 = @(x) sign(x)./( 1+abs(x) ).^2;
                    % g3' = -2/(1+abs(x))^3
                    grad = X * (sign_y .* inv_pabs_y .* inv_pabs_y)' / T; %grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
                    if fastica_on
                        Edgx = sum(-2*inv_pabs_y.*inv_pabs_y.*inv_pabs_y) / T; % Edgx = sum(-2*inv_pabs_y.^3)/T;
                    else
                        vEGx = 2 * (EGx(3) > nf3.critical_point) - 1;
                    end

                case 5

                    % G5 = @(x) x*|x|/(10+|x|);
                    % g5 = @(x) |x|*(20+|x|)/(10+|x|)^2;
                    % g5 = 200*sign(x) / (10+|x|)^3
                    grad = X * (abs_y .* (20 + abs_y) .* inv_p10abs_y.^2)' / T;
                    if fastica_on
                        Edgx = sum(200*sign_y.*inv_p10abs_y.*inv_p10abs_y.*inv_p10abs_y) / T; %Edgx = sum( 200*sign_y.*inv_p10abs_y.^3 )/T;
                    else
                        vEGx = 2 * (EGx(5) > nf5.critical_point) - 1;
                    end

                case 7

                    % G7 = @(x) x/(1+x^2)
                    % g7 = @(x) (1-x^2)/(1+x^2)^2
                    % g7' = 2x(x^2-3)/(1+x^2)^3
                    grad = X * ((1 - yy) .* inv_pabs_yy.^2)' / T;
                    if fastica_on
                        Edgx = sum(2*y.*(yy - 3).*inv_pabs_yy.*inv_pabs_yy.*inv_pabs_yy) / T; %Edgx = sum( 2*y.*(yy-3).*inv_pabs_yy.^3 )/T;
                    else
                        vEGx = 2 * (EGx(7) > nf7.critical_point) - 1;
                    end

                otherwise
                    ;
            end

            if fastica_on
                w1 = grad - Edgx * w;
            else
                grad = vEGx * grad;
                grad = grad - (w' * grad) * w;
                grad = grad / norm(grad);
                w1 = w + mu * grad;
            end

            W(n, :) = w1';

        end

        W = symdecor(W);

        if Cost(iter) < min_cost
            cost_increase_counter = 0;
            min_cost = Cost(iter);
            best_W = last_W;
            max_negentropy = negentropy_array;
        else
            cost_increase_counter = cost_increase_counter + 1;
        end
        min_cost_queue(iter) = min_cost;

        if fastica_on
            if cost_increase_counter >= max_cost_increase_number || 1 - min(abs(diag(W*last_W'))) < tolerance
                cost_increase_counter = 0;
                W = best_W;
                last_W = W;
                iter_fastica = iter;
                fastica_on = 0;
                continue;
            end
        else
            if cost_increase_counter > max_cost_increase_number
                break;
            end
        end

        last_W = W;

    end
    W = best_W;
    if show_cost
        subplot(1, 4, 3);
        plot([1:iter_fastica], Cost([1:iter_fastica]), 'r');
        hold on;
        plot([iter_fastica + 1:iter], Cost([iter_fastica + 1:iter]), 'b');
        xlabel('Number of iterations')
        ylabel('Cost')
        title('(b) Refine orthogonal ICA')
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    if show_cost
        subplot(1, 4, 3);
        title('(b) No saddles, no Refinement')
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort all the component
% [max_negentropy, index_sort] = sort(max_negentropy, 'descend');
% W = W(index_sort, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Part IV: non-orth ICA
%   fixed small step size for refinement, gradient search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    fprintf('\nNonorthogonal ICA stage.');
end

% have to realign W since there is a sorting of W based on NG in earlier
% section
num_guess = size(guess_mat, 2);
num_W = size(W, 1);
%Re-sort existing W based on correlation with matrix W:
for kl = 1:num_guess
    r_n_c = guess_mat(:, kl);
    for lp = 1:num_W
        w = W(lp, :).';
        R_par = corrcoef(X'*w, r_n_c); %corr(y',rnc)
        corr_w_guess(kl, lp) = R_par(1, 2);
    end
end
%We may need to use auction to chose order:
[~, max_idx] = max(abs(corr_w_guess), [], 2); % ; Added by Zois

if length(unique(max_idx)) ~= num_guess
    [colsol, ~] = auction((1 - abs(corr_w_guess))');
    max_idx = colsol';
end

c = setxor((1:num_W).', max_idx);
sort_order = [max_idx; c];

W = W(sort_order, :);


last_W = W;
best_W = W;
Cost = zeros(max_iter_nonorth, 1);
% min_cost = inf;
min_cost_queue = min_cost * ones(max_iter_nonorth, 1);
error = inf;
mu = 1; %step size
min_mu = 1e-6;
max_cost_increase_number = 3;
cost_increase_counter = 0;
mu_idx = false(size(mu_c));

decaying_factor = 0.95;
min_change = 1e-6;
min_iter = 100;
mu_old = mu_c;

for iter = 1:max_iter_nonorth

    Cost(iter) = -log(abs(det(W)));
    %W = W(randperm(N),:);
    for n = 1:N

        % if N is too large, we adopt a smarter way to calculate h_n to
        % reduce the computational load
        if N > 7

            if n == 1
                Wn = W(2:N, :);
                inv_Q = inv(Wn*Wn');
            else
                n_last = n - 1;
                Wn_last = [W(1:n_last-1, :); W(n_last+1:N, :)];
                w_current = W(n, :)';
                w_last = W(n_last, :)';
                c = Wn_last * (w_last - w_current);
                c(n_last) = 0.5 * (w_last' * w_last - w_current' * w_current);
                e_last = zeros(N-1, 1);
                e_last(n_last) = 1;

                temp1 = inv_Q * c;
                temp2 = inv_Q(:, n_last);
                inv_Q_plus = inv_Q - temp1 * temp2' / (1 + temp1(n_last));

                temp1 = inv_Q_plus' * c;
                temp2 = inv_Q_plus(:, n_last);
                inv_Q = inv_Q_plus - temp2 * temp1' / (1 + c' * temp2);
                % inv_Q is Hermitian
                inv_Q = (inv_Q + inv_Q') / 2;
            end

            temp1 = randn(N, 1);
            W_n = [W(1:n-1, :); W(n+1:N, :)];
            h = temp1 - W_n' * inv_Q * W_n * temp1;

        else
            temp1 = randn(N, 1);
            temp2 = [W(1:n-1, :); W(n+1:N, :)];
            h = temp1 - temp2' * inv(temp2*temp2') * temp2 * temp1;
        end

        w = W(n, :)';
        y = w' * X;

        % evaluate the upper bound of negentropy of the n-th component
        NE_Bound = zeros(K, 1);
        EGx = zeros(K, 1);

        % we only need to calculate these quantities once
        yy = y .* y;
        sign_y = sign(y);
        abs_y = sign_y .* y;
        inv_pabs_y = 1 ./ (1 + abs_y);
        inv_pabs_yy = 1 ./ (1 + yy);
        inv_p10abs_y = 1 ./ (10 + abs_y);

        % G1(x) = x^4
        EGx(1) = sum(yy.*yy) / T;
        if EGx(1) < nf1.min_EGx
            NE_Bound(1) = simplified_ppval(nf1.pp_slope, nf1.min_EGx) * (EGx(1) - nf1.min_EGx);
            NE_Bound(1) = simplified_ppval(nf1.pp, nf1.min_EGx) + abs(NE_Bound(1));
        else if EGx(1) > nf1.max_EGx
                % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                NE_Bound(1) = 0;
        else
            NE_Bound(1) = simplified_ppval(nf1.pp, EGx(1));
        end
        end

        % G3 = @(x) abs(x)./( 1+abs(x) );
        EGx(3) = 1 - sum(inv_pabs_y) / T;
        if EGx(3) < nf3.min_EGx
            NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.min_EGx) * (EGx(3) - nf3.min_EGx);
            NE_Bound(3) = simplified_ppval(nf3.pp, nf3.min_EGx) + abs(NE_Bound(3));
        else if EGx(3) > nf3.max_EGx
                NE_Bound(3) = simplified_ppval(nf3.pp_slope, nf3.max_EGx) * (EGx(3) - nf3.max_EGx);
                NE_Bound(3) = simplified_ppval(nf3.pp, nf3.max_EGx) + abs(NE_Bound(3));
        else
            NE_Bound(3) = simplified_ppval(nf3.pp, EGx(3));
        end
        end

        % G5 = @(x) x*|x|/(10+|x|);
        EGx(5) = sum(y.*abs_y.*inv_p10abs_y) / T;
        if EGx(5) < nf5.min_EGx
            NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.min_EGx) * (EGx(5) - nf5.min_EGx);
            NE_Bound(5) = simplified_ppval(nf5.pp, nf5.min_EGx) + abs(NE_Bound(5));
        else if EGx(5) > nf5.max_EGx
                NE_Bound(5) = simplified_ppval(nf5.pp_slope, nf5.max_EGx) * (EGx(5) - nf5.max_EGx);
                NE_Bound(5) = simplified_ppval(nf5.pp, nf5.max_EGx) + abs(NE_Bound(5));
        else
            NE_Bound(5) = simplified_ppval(nf5.pp, EGx(5));
        end
        end

        % G7 = @(x) x/(1+x^2);
        EGx(7) = sum(y.*inv_pabs_yy) / T;
        if EGx(7) < nf7.min_EGx
            NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.min_EGx) * (EGx(7) - nf7.min_EGx);
            NE_Bound(7) = simplified_ppval(nf7.pp, nf7.min_EGx) + abs(NE_Bound(7));
        else if EGx(7) > nf7.max_EGx
                NE_Bound(7) = simplified_ppval(nf7.pp_slope, nf7.max_EGx) * (EGx(7) - nf7.max_EGx);
                NE_Bound(7) = simplified_ppval(nf7.pp, nf7.max_EGx) + abs(NE_Bound(7));
        else
            NE_Bound(7) = simplified_ppval(nf7.pp, EGx(7));
        end
        end

        [max_NE, max_i] = max(NE_Bound); % select the tightest bound
        Cost(iter) = Cost(iter) - max_NE;
        %

        %         alpha = randn(n_refs, K);
        %         mu = max(0, alpha);
        %         mu_idx = false(size(mu));
        %         idx = arrayfun(@(e)(find(Params.rhoList > e, 1)), abs(eps_yr)); % MxK
        %         idx(mu_idx) = idx(mu_idx) - 1;
        %         rho = Params.rhoList(idx); % MxK
        %         mu = max(0, alpha);
        %         mu_idx = mu_idx | (mu >= Params.maxMu);
        %         mu_idx = mu_idx & (mu > 0);
        %         mu = min(Params.maxMu, mu);
        %         alpha = Params.gamma * (rho - abs(eps_yr)) + mu;

        %

        % CONSTRAINT
        if n <= num_guess
            %constrain reference signals
            r_n_c = guess_mat(:, n);
            e_pair = corrcoef(y', r_n_c); %correlation of y and r_n
            %e_pair_s = corrcoef(S(n,:)',r_n_c); %correlation of x_n and r_n
            dis_wr = abs(e_pair(1, 2));
            %dis_wr_true = abs(e_pair_s(1,2))
            %dis_wr_true_arr(iter,n) = dis_wr_true;
            if mu_idx(n)
                rho_n = max(rho(rho < dis_wr));
            else
                rho_n = min(rho(rho > dis_wr));
            end
            if isempty(rho_n)
                rho_n = 0.01;
            end

            rho_n_arr(iter, n) = rho_n;
            %            [min_val,ind_min] = min(abs(rho-dis_wr));
            %            rho_n = rho(ind_min)

            mu_old(n) = mu_c(n);
            mu_idx(n) = mu_idx(n) | (mu_c(n) >= 1);
            mu_idx(n) = mu_idx(n) & (mu_c(n) > 0);
            mu_c(n) = min(1, mu_c(n));
            mu_c(n) = max(0, mu_c(n)+gam*(rho_n - dis_wr));
            mu_signed(n) = sign(e_pair(1, 2)) * mu_c(n);
            
            r_n_c = r_n_c / norm(r_n_c);
            dis_wr_all(iter, n) = dis_wr;
        end

        switch max_i

            case 1

                %                     % G1(x) = x^4
                %                     % g1(x) = 4x^3
                %                     % g1' =  12x^2
                %                     vEGx = 2*( EGx(1) > nf1.critical_point ) - 1;
                %                     grad = X*( 4*y.*yy )'/T;
                EGx(1) = max(min(EGx(1), nf1.max_EGx), nf1.min_EGx);
                grad = h / (h' * w) + X * (4 * y .* yy)' * simplified_ppval(nf1.pp_slope, EGx(1)) / T; % gradient

            case 3

                %                     % G3 = @(x) abs(x)./( 1+abs(x) );
                %                     % g3 = @(x) sign(x)./( 1+abs(x) ).^2;
                %                     % g3' = -2/(1+abs(x))^3
                %                     vEGx = 2*( EGx(3) > nf3.critical_point ) - 1;
                %                     grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
                EGx(3) = max(min(EGx(3), nf3.max_EGx), nf3.min_EGx);
                grad = h / (h' * w) + X * (sign_y .* inv_pabs_y.^2)' * simplified_ppval(nf3.pp_slope, EGx(3)) / T;

            case 5

                %                     % G5 = @(x) x*|x|/(10+|x|);
                %                     % g5 = @(x) |x|*(20+|x|)/(10+|x|)^2;
                %                     % g5 = 200*sign(x) / (10+|x|)^3
                %                     vEGx = 2*( EGx(5) > nf5.critical_point ) - 1;
                %                     grad = X*( abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'/T;
                EGx(5) = max(min(EGx(5), nf5.max_EGx), nf5.min_EGx);
                grad = h / (h' * w) + X * (abs_y .* (20 + abs_y) .* inv_p10abs_y.^2)' * simplified_ppval(nf5.pp_slope, EGx(5)) / T;

            case 7

                %                     % G7 = @(x) x/(1+x^2)
                %                     % g7 = @(x) (1-x^2)/(1+x^2)^2
                %                     vEGx = 2*( EGx(7) > nf7.critical_point ) - 1;
                %                     grad = X*( (1-yy).*inv_pabs_yy.^2 )'/T;
                EGx(7) = max(min(EGx(7), nf7.max_EGx), nf7.min_EGx);
                grad = h / (h' * w) + X * ((1 - yy) .* inv_pabs_yy.^2)' * simplified_ppval(nf7.pp_slope, EGx(7)) / T; %this is the negative value of eq14

            otherwise
                ;
        end
        if n <= num_guess
            grad = grad + mu_signed(n) .* (X * r_n_c) * (1 / sqrt(T)); % A column vector. Gradient should be the original one and the constraint one
        end
        grad = grad - (w' * grad) * w;
        grad = grad / norm(grad);

        w1 = w + mu * grad;
        w1 = w1 / norm(w1);
        W(n, :) = w1';


    end
    %     Cost(iter) = Cost(iter) - sum((mu_new.^2 - mu_old.^2))/(2*gam);
    %
    %     if Cost(iter) < min_cost
    %         cost_increase_counter = 0;
    %         min_cost = Cost(iter);
    %         best_W = last_W;
    %     else
    %         cost_increase_counter = cost_increase_counter + 1;
    %     end
    %     min_cost_queue(iter) = min_cost;
    %
    %     if cost_increase_counter > max_cost_increase_number
    %         if mu > min_mu
    %             cost_increase_counter = 0;
    %             W = best_W;
    %             last_W = W;
    %             mu = mu/2;
    %             continue;
    %         else
    %             currentChange = max(0,max(1-abs(diag(last_W*W'))))
    %             if currentChange < 1e-6
    %                 break;
    %             end
    %         end
    %     else
    %         last_W = W;
    %     end

    Cost(iter) = Cost(iter) - sum((mu_c.^2 - mu_old.^2)) / (2 * gam);
    %mu = max(mu * decaying_factor,min_mu)
    mu = max((decaying_factor^iter)*1, min_mu);
    currentChange = max(0, max(1-abs(diag(last_W*W'))));
    if currentChange < min_change && iter > min_iter
        best_W = W;
        break;
    else
        last_W = W;
    end

end
W = best_W;
W = W * P;

% Cyrus 060424 changing signs of components for optimization below
S = W*xOrig;
for i=1:length(W)
    nixIcSignFlip1(i)=corr(S(i,:)',guess_mat(:,i));
end
nixIcSignFlip2 = find(nixIcSignFlip1 < 0);
W(nixIcSignFlip2,:) = W(nixIcSignFlip2,:)*-1; % Flipping on W according with matrix algebra
% Cyrus 060424 changing signs of components for optimization above

if show_cost
    subplot(1, 4, 4);
    plot([1:iter], Cost([1:iter]));
    xlabel('Number of iterations')
    ylabel('Cost')
    title('(c) Nonorthogonal ICA')
end
% [EOF ICA]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this subfunction is faster than ppval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = simplified_ppval(pp, xs)
% a simplified version of ppval

b = pp.breaks;
c = pp.coefs;
l = pp.pieces;
k = 4; % k = pp.order;
dd = 1; % dd = pp.dim;

% find index
index = 0;
middle_index = 0;
if xs > b(l)
    index = l;
else if xs < b(2)
        index = 1;
else

    low_index = 1;
    high_index = l;
    while 1

        middle_index = round(0.6*low_index+0.4*high_index);
        if b(middle_index) > xs
            high_index = middle_index;
        else
            low_index = middle_index;
        end

        if low_index == high_index - 1
            index = low_index;
            break;
        end

    end

end
end


% now go to local coordinates ...
xs = xs - b(index);

% ... and apply nested multiplication:
v = c(index, 1);
for i = 2:k
    v = xs * v + c(index, i);
end


function W = symdecor(M)
%fast symmetric orthogonalization
[V, D] = eig(M*M');
W = (V .* (ones(size(M, 2), 1) * (1 ./ sqrt(diag(D)')))) * V' * M;


function [X, P] = pre_processing(X)
% pre-processing program
[N, T] = size(X);
% remove DC
Xmean = mean(X, 2);
X = X - Xmean * ones(1, T);

% spatio pre-whitening 1
R = X * X' / T;
P = inv_sqrtmH(R); %P = inv(sqrtm(R));
X = P * X;


function A = inv_sqrtmH(B)
%
[V, D] = eig(B);
d = diag(D);
d = 1 ./ sqrt(d);
A = V * diag(d) * V';

function properties = getopt(properties, varargin)
%GETOPT - Process paired optional arguments as 'prop1',val1,'prop2',val2,...
%
%   getopt(properties,varargin) returns a modified properties structure,
%   given an initial properties structure, and a list of paired arguments.
%   Each argumnet pair should be of the form property_name,val where
%   property_name is the name of one of the field in properties, and val is
%   the value to be assigned to that structure field.
%
%   No validation of the values is performed.

%%
% EXAMPLE:
%   properties = struct('zoom',1.0,'aspect',1.0,'gamma',1.0,'file',[],'bg',[]);
%   properties = getopt(properties,'aspect',0.76,'file','mydata.dat')
% would return:
%   properties =
%         zoom: 1
%       aspect: 0.7600
%        gamma: 1
%         file: 'mydata.dat'
%           bg: []
%
% Typical usage in a function:
%   properties = getopt(properties,varargin{:})

% Function from
% http://mathforum.org/epigone/comp.soft-sys.matlab/sloasmirsmon/bp0ndp$crq5@cui1.lmms.lmco.com

% dgleich
% 2003-11-19
% Added ability to pass a cell array of properties

if ~isempty(varargin) && (iscell(varargin{1}))
    varargin = varargin{1};
end

% Process the properties (optional input arguments)
prop_names = fieldnames(properties);
TargetField = [];
for ii = 1:length(varargin)
    arg = varargin{ii};
    if isempty(TargetField)
        if ~ischar(arg)
            error('Property names must be character strings');
        end
        %f = find(strcmp(prop_names, arg));
        if isempty(find(strcmp(prop_names, arg), 1)) %length(f) == 0
            error('%s ', ['invalid property ''', arg, '''; must be one of:'], prop_names{:});
        end
        TargetField = arg;
    else
        properties.(TargetField) = arg;
        TargetField = '';
    end
end
if ~isempty(TargetField)
    error('Property names and values must be specified in pairs.');
end

function [colsol, rowsol] = auction(assignCost)

%function [colsol, rowsol] = auction(assignCost, guard)
%AUCTION : performs assignment using Bertsekas' auction algorithm
%          The auction is 1-sided and uses a fixed epsilon.
%  colsol = auction(assignCost) returns column assignments: colsol(j) gives the
%           row assigned to column j.
%
%           assignCost is an m x n matrix of costs for associating each row
%           with each column.  m >= n (rectangular assignment is allowed).
%           Auction finds the assignment that minimizes the costs.
%
%  colsol = auction(assignCost, guard) sets the cost of column
%           non-assignment = guard.  All assignments will have cost < guard.
%           A column will not be assigned if the savings from eliminating it
%           exceeds the guard cost.  colsol(j) = 0 indicates that the optimal
%           solution does not assign column j to any row.
%
%  [colsol,rowsol] = auction(assignCost) also returns the row assignments:
%           rowsol(i) gives the column assigned to row j.
%
%  Reference
%  Bertsekas, D. P., "The Auction Algorithm: A Distributed Relaxation Method
%  for the Assignment Problem," Annals of Operations Research, 14, 1988, pp.
%  103-123.

% Mark Levedahl
%
%=================================================================================

[m, n] = size(assignCost);

if m < n
    error('cost matrix must have no more columns than rows.')
end

% augment cost matrix with a guard row if specified.
m0 = m;
% if isfinite(guard) %%% DONE BY JS
% 	m = m+1;
% 	assignCost(m,:) = guard;
% end

% init return arrays
colsol = zeros(1, n);
rowsol = zeros(m, 1);
price = zeros(m, 1);
EPS = sqrt(eps) / (n + 1);

% 1st step is a full parallel solution.  Get bids for all columns
jp = 1:n;
f = assignCost;
[b1, ip] = min(f); % cost and row of best choice
f(ip+m*(0:n - 1)) = inf; % eliminate best from contention
bids = min(f) - b1; % cost of runner up choice hence bid

% Arrange bids so highest (best) are last and will overwrite the lower bids.
[tmp, ibid] = sort(bids(:));

% Now assign best bids (lesser bids are overwritten by better ones).
price(ip(ibid)) = price(ip(ibid)) + EPS + tmp;
rowsol(ip(ibid)) = jp(ibid); % do row assignments
iy = find(rowsol);
colsol(rowsol(iy)) = iy; % map to column assignments

% The guard row cannot be assigned (always available)
if m0 < m
    price(m) = 0;
    rowsol(m) = 0;
end

% Now Continue with non-parallel code handling any contentions.
while ~all(colsol)
    for jp = find(~colsol)
        f = assignCost(:, jp) + price; % costs
        [b1, ip] = min(f); % cost and row of best choice
        if ip > m0
            colsol(jp) = m;
        else
            f(ip) = inf; % eliminate from contention
            price(ip) = price(ip) + EPS + min(f) - b1; % runner up choice hence bid
            if rowsol(ip) % take the row away if already assigned
                colsol(rowsol(ip)) = 0;
            end
            rowsol(ip) = jp; % update row and column assignments
            colsol(jp) = ip;
        end % if ip == m
    end
end

% screen out infeasible assignments
if m > m0
    colsol(colsol == m) = 0;
    rowsol(m) = [];
end
return
