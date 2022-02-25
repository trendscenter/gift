function W = icatb_constrained_ICA_EBM(X,guess_mat,rho,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% Program by Xi-Lin Li. Please contact me at lixilin@umbc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% References:
% 
% [1] Xi-Lin Li and Tulay Adali, "A novel entropy estimator and its application to ICA," 
% IEEE International Workshop on Machine Learning for Signal Processing 2009. 
%
% [2] Xi-Lin Li and Tulay Adali, "Independent component analysis by entropy bound minimization"
% IEEE Transaction on Signal Processing, submitted.
% 
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Part 0: Here begins some pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% build default options structure
options=struct( ...
   'opt_approach','newton', ... % optimization type: gradient, (newton), quasi
   'complex_valued',false, ... % if any input is complex or quasi approach used then setting is forced to true
   'circular',false, ... % set to true to only consider circular for complex-valued cost function
   'whiten',true, ... % whitening is optional (except for quasi approach it is required)
   'verbose',false, ... % verbose true enables print statements
   'A',[], ... % true mixing matrices A, automatically sets verbose
   'W_init',[], ... % initial estimates for demixing matrices in W
   'jdiag_initW',false, ... % use CCA (K=2) or joint diagonalization (K>2)
   'maxIter',512, ... % max number of iterations
   'WDiffStop',1e-6, ... % stopping criterion
   'alpha0',1.0, ... % initial step size scaling (will be doubled for complex-valued)
   'Save_W',false ...
   );

% load in user supplied options
options=getopt(options,varargin{:});

supplyA=~isempty(options.A); % set to true if user has supplied true mixing matrices
blowup = 1e3;
alphaScale=0.9; % alpha0 to alpha0*alphaScale when cost does not decrease
alphaMin=options.WDiffStop; % alpha0 will max(alphaMin,alphaScale*alpha0)
outputISI=false;
opt_approach=find(strcmp(options.opt_approach, ...
   {'gradient','newton','quasi'}),1); % 1='gradient', 2='newton', 3='quasi'
% options.whiten=(options.whiten || (opt_approach==3)); % whitening required for real-valued quasi & complex-valued gradient approach
% filename = 'IVAG_results.txt';



max_iter_fastica = 100;
max_iter_orth = 1000;
max_iter_orth_refine = 1000;
max_iter_nonorth = 1000;
saddle_test_enable = 1;
tolerance = 1e-4;
max_cost_increase_number = 5;
stochastic_search_factor = 1;
% mu_c = 0.5;
mu_c = 0;

gam = 3;
% rho = rho;
num_guess=size(guess_mat,2);

verbose = 0;        % report the progress if verbose==1
show_cost = 0;      % show the cost values vs. iterations at each stage if show_cost==1

% Load 8 measuring functions. But we only use 4 of them.
K = 8;          
load icatb_nf_table; 

[N,T] = size(X);
[X, P] = pre_processing( X );
%% Initialize W
if ~isempty(options.W_init)
   W=options.W_init;
else
   
   % use SEA, or fastica with pow3, to provide an initial guess
    W = randn(N,N);
end

W = symdecor(W);
last_W = W;
best_W = W;
Cost = zeros(max_iter_fastica,1);
min_cost = inf;
cost_increase_counter = 0;
for iter = 1 : max_iter_fastica
    
    Y = W*X;
    
    for n = 1 : N
        
        y = Y(n,:);     %y = w'*X;
        
        % evaluate the upper bound of negentropy of the nth component
        NE_Bound = zeros(K,1);
        EGx = zeros(K,1);
        
        % we only need to calculate these quantities once
        yy = y.*y;                              
        sign_y = sign(y);
        abs_y = sign_y.*y;
        inv_pabs_y = 1./(1+abs_y);
        inv_pabs_yy = 1./(1+yy);
        inv_p10abs_y = 1./(10+abs_y);
              
        % G1(x) = x^4
        EGx(1) = sum( yy.*yy )/T;  
        if EGx(1)<nf1.min_EGx
            NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.min_EGx ) * ( EGx(1) - nf1.min_EGx );
            NE_Bound(1) = simplified_ppval( nf1.pp, nf1.min_EGx ) + abs( NE_Bound(1) );
        else if EGx(1)>nf1.max_EGx
                % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                NE_Bound(1) = 0;
            else
                NE_Bound(1) =  simplified_ppval( nf1.pp, EGx(1) );
            end
        end
        
        % G3 = @(x) abs(x)./( 1+abs(x) );
        EGx(3) = 1-sum( inv_pabs_y )/T;
        if EGx(3)<nf3.min_EGx
            NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.min_EGx ) * ( EGx(3) - nf3.min_EGx );
            NE_Bound(3) = simplified_ppval( nf3.pp, nf3.min_EGx ) + abs( NE_Bound(3) );
        else if EGx(3)>nf3.max_EGx
                NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.max_EGx ) * ( EGx(3) - nf3.max_EGx );
                NE_Bound(3) = simplified_ppval( nf3.pp, nf3.max_EGx ) + abs( NE_Bound(3) );
            else
                NE_Bound(3) =  simplified_ppval( nf3.pp, EGx(3) );
            end
        end
        
        % G5 = @(x) x*|x|/(10+|x|);
        EGx(5) = sum( y.*abs_y.*inv_p10abs_y )/T;   
        if EGx(5)<nf5.min_EGx
            NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.min_EGx ) * ( EGx(5) - nf5.min_EGx );
            NE_Bound(5) = simplified_ppval( nf5.pp, nf5.min_EGx ) + abs( NE_Bound(5) );
        else if EGx(5)>nf5.max_EGx
                NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.max_EGx ) * ( EGx(5) - nf5.max_EGx );
                NE_Bound(5) = simplified_ppval( nf5.pp, nf5.max_EGx ) + abs( NE_Bound(5) );
            else
                NE_Bound(5) =  simplified_ppval( nf5.pp, EGx(5) );
            end
        end
        
        % G7 = @(x) x/(1+x^2);
        EGx(7) = sum( y.*inv_pabs_yy )/T;   
        if EGx(7)<nf7.min_EGx
            NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.min_EGx ) * ( EGx(7) - nf7.min_EGx );
            NE_Bound(7) = simplified_ppval( nf7.pp, nf7.min_EGx ) + abs( NE_Bound(7) );
        else if EGx(7)>nf7.max_EGx
                NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.max_EGx ) * ( EGx(7) - nf7.max_EGx );
                NE_Bound(7) = simplified_ppval( nf7.pp, nf7.max_EGx ) + abs( NE_Bound(7) );
            else
                NE_Bound(7) =  simplified_ppval( nf7.pp, EGx(7) );
            end
        end
               
        % select the tightest bound
        [max_NE, max_i] = max( NE_Bound );  
        negentropy_array(n) = max_NE;
        Cost(iter) = Cost(iter) - max_NE;
        
    end
    
    if Cost(iter)<min_cost
        min_cost = Cost(iter);
        best_W = last_W;
        cost_increase_counter = 0;
    else
        cost_increase_counter = cost_increase_counter + 1;
    end
    
    W = (Y.*Y.*Y)*X'/T - 3*W;
    W = symdecor(W);
    
    if 1-min(abs(diag(W*last_W'))) < tolerance
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
    figure; subplot(1,4,1); plot(1:iter,Cost(1:iter));
    xlabel('Number of iterations')
    ylabel('Cost')
    title('(0) Initialize use SEA')
end

% [End of Part 0]
i=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Part I: orthogonal ICA
%  varying step size, stochastic gradient search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    fprintf('\nOrthogonal ICA stage.');
end
last_W = W;
best_W = W;
Cost = zeros(max_iter_orth,1);
min_cost = inf;
min_cost_queue = min_cost*ones(max_iter_orth,1);
mu = 1/6.25;
min_mu = 1/50;
cost_increase_counter = 0;
fastica_on = 1;
error = 0;
max_negentropy = zeros(N,1);
negentropy_array = zeros(N,1);
for iter = 1 : max_iter_orth
    
    Y = W*X;
    for n = 1 : N
        
        w = W(n,:)';
        y = Y(n,:);     %y = w'*X;
        
        % evaluate the upper bound of negentropy of the nth component
        NE_Bound = zeros(K,1);
        EGx = zeros(K,1);
        
        % we only need to calculate these quantities once
        yy = y.*y;                              
        sign_y = sign(y);
        abs_y = sign_y.*y;
        inv_pabs_y = 1./(1+abs_y);
        inv_pabs_yy = 1./(1+yy);
        inv_p10abs_y = 1./(10+abs_y);
              
        % G1(x) = x^4
        EGx(1) = sum( yy.*yy )/T;  
        if EGx(1)<nf1.min_EGx
            NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.min_EGx ) * ( EGx(1) - nf1.min_EGx );
            NE_Bound(1) = simplified_ppval( nf1.pp, nf1.min_EGx ) + abs( NE_Bound(1) );
        else if EGx(1)>nf1.max_EGx
                % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                NE_Bound(1) = 0;
            else
                NE_Bound(1) =  simplified_ppval( nf1.pp, EGx(1) );
            end
        end
        
        % G3 = @(x) abs(x)./( 1+abs(x) );
        EGx(3) = 1-sum( inv_pabs_y )/T;
        if EGx(3)<nf3.min_EGx
            NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.min_EGx ) * ( EGx(3) - nf3.min_EGx );
            NE_Bound(3) = simplified_ppval( nf3.pp, nf3.min_EGx ) + abs( NE_Bound(3) );
        else if EGx(3)>nf3.max_EGx
                NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.max_EGx ) * ( EGx(3) - nf3.max_EGx );
                NE_Bound(3) = simplified_ppval( nf3.pp, nf3.max_EGx ) + abs( NE_Bound(3) );
            else
                NE_Bound(3) =  simplified_ppval( nf3.pp, EGx(3) );
            end
        end
        
        % G5 = @(x) x*|x|/(10+|x|);
        EGx(5) = sum( y.*abs_y.*inv_p10abs_y )/T;   
        if EGx(5)<nf5.min_EGx
            NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.min_EGx ) * ( EGx(5) - nf5.min_EGx );
            NE_Bound(5) = simplified_ppval( nf5.pp, nf5.min_EGx ) + abs( NE_Bound(5) );
        else if EGx(5)>nf5.max_EGx
                NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.max_EGx ) * ( EGx(5) - nf5.max_EGx );
                NE_Bound(5) = simplified_ppval( nf5.pp, nf5.max_EGx ) + abs( NE_Bound(5) );
            else
                NE_Bound(5) =  simplified_ppval( nf5.pp, EGx(5) );
            end
        end
        
        % G7 = @(x) x/(1+x^2);
        EGx(7) = sum( y.*inv_pabs_yy )/T;   
        if EGx(7)<nf7.min_EGx
            NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.min_EGx ) * ( EGx(7) - nf7.min_EGx );
            NE_Bound(7) = simplified_ppval( nf7.pp, nf7.min_EGx ) + abs( NE_Bound(7) );
        else if EGx(7)>nf7.max_EGx
                NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.max_EGx ) * ( EGx(7) - nf7.max_EGx );
                NE_Bound(7) = simplified_ppval( nf7.pp, nf7.max_EGx ) + abs( NE_Bound(7) );
            else
                NE_Bound(7) =  simplified_ppval( nf7.pp, EGx(7) );
            end
        end
               
        % select the tightest bound
        [max_NE, max_i] = max( NE_Bound );  
        negentropy_array(n) = max_NE;
        Cost(iter) = Cost(iter) - max_NE;
        
        if ~fastica_on
            weight = rand(1,T);
        end
                   
            % Perform orthogonal ICA
            switch max_i
                
                case 1
                    
                    % G1(x) = x^4
                    % g1(x) = 4x^3
                    % g1' =  12x^2
                    if fastica_on 
                        grad = X*( 4*y.*yy )'/T;
                        Edgx = 12; 
                    else
                        grad = X*( 4*weight.*y.*yy )'/sum(weight);
                        vEGx = 2*( EGx(1) > nf1.critical_point ) - 1;
                    end
                    
                case 3
                    
                    % G3 = @(x) abs(x)./( 1+abs(x) );
                    % g3 = @(x) sign(x)./( 1+abs(x) ).^2;
                    % g3' = -2/(1+abs(x))^3
                    if fastica_on 
                        grad = X*( sign_y.*inv_pabs_y.*inv_pabs_y )'/T;      %grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
                        Edgx = sum(-2*inv_pabs_y.*inv_pabs_y.*inv_pabs_y)/T; % Edgx = sum(-2*inv_pabs_y.^3)/T; 
                    else
                        grad = X*( weight.*sign_y.*inv_pabs_y.*inv_pabs_y )'/sum(weight);      %grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
                        vEGx = 2*( EGx(3) > nf3.critical_point ) - 1;
                    end
                    
                case 5
                    
                    % G5 = @(x) x*|x|/(10+|x|);
                    % g5 = @(x) |x|*(20+|x|)/(10+|x|)^2;      
                    % g5 = 200*sign(x) / (10+|x|)^3
                    if fastica_on 
                        grad = X*( abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'/T;
                        Edgx = sum( 200*sign_y.*inv_p10abs_y.*inv_p10abs_y.*inv_p10abs_y )/T;    %Edgx = sum( 200*sign_y.*inv_p10abs_y.^3 )/T; 
                    else
                        grad = X*( weight.*abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'/sum(weight);
                        vEGx = 2*( EGx(5) > nf5.critical_point ) - 1;
                    end
                    
                case 7
                    
                    % G7 = @(x) x/(1+x^2)
                    % g7 = @(x) (1-x^2)/(1+x^2)^2
                    % g7' = 2x(x^2-3)/(1+x^2)^3
                    if fastica_on 
                        grad = X*( (1-yy).*inv_pabs_yy.^2 )'/T;
                        Edgx = sum( 2*y.*(yy-3).*inv_pabs_yy.*inv_pabs_yy.*inv_pabs_yy )/T;    %Edgx = sum( 2*y.*(yy-3).*inv_pabs_yy.^3 )/T; 
                    else
                        grad = X*( weight.*(1-yy).*inv_pabs_yy.^2 )'/sum(weight);
                        vEGx = 2*( EGx(7) > nf7.critical_point ) - 1;
                    end
                                        
                otherwise
                    ;
            end
            
            if fastica_on
                w1 = grad - Edgx*w;
            else  
                grad = vEGx*grad;
                grad = grad - (w'*grad)*w;
                grad = grad / norm(grad);
                w1 = w + mu*grad;
            end

            W(n,:) = w1';
 
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
        if cost_increase_counter>=max_cost_increase_number || 1-min(abs(diag(W*last_W'))) < tolerance
            cost_increase_counter = 0;
            W = best_W;
            last_W = W;
            iter_fastica = iter;
            fastica_on = 0;
            continue;
        end
    else
        if cost_increase_counter > stochastic_search_factor*max_cost_increase_number
            if mu > min_mu
                cost_increase_counter = 0;
                W = best_W;
                last_W = W;
                mu = mu/2;
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
    subplot(1,4,2);
    plot([1:iter_fastica], Cost([1:iter_fastica]),'r');
    hold on; plot([iter_fastica+1:iter], Cost([iter_fastica+1:iter]),'b');
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
    Y = W*X;
    for m = 1 : N
        w1 = W(m,:)';
        ym = Y(m,:);    %ym = w1'*X;
        for n = m+1 : N
            w2 = W(n,:)';
            yn = Y(n,:);    %yn = w2'*X;
            
            yr1 = (ym+yn)/sqrt(2);
            yr2 = (ym-yn)/sqrt(2);
            
            y = yr1;
            % we only need to calculate these quantities once
            yy = y.*y;                              
            sign_y = sign(y);
            abs_y = sign_y.*y;
            inv_pabs_y = 1./(1+abs_y);
            inv_pabs_yy = 1./(1+yy);
            inv_p10abs_y = 1./(10+abs_y);

            % G1(x) = x^4
            EGx(1) = sum( yy.*yy )/T;  
            if EGx(1)<nf1.min_EGx
                NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.min_EGx ) * ( EGx(1) - nf1.min_EGx );
                NE_Bound(1) = simplified_ppval( nf1.pp, nf1.min_EGx ) + abs( NE_Bound(1) );
            else if EGx(1)>nf1.max_EGx
                    % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                    % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                    NE_Bound(1) = 0;
                else
                    NE_Bound(1) =  simplified_ppval( nf1.pp, EGx(1) );
                end
            end

            % G3 = @(x) abs(x)./( 1+abs(x) );
            EGx(3) = 1-sum( inv_pabs_y )/T;
            if EGx(3)<nf3.min_EGx
                NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.min_EGx ) * ( EGx(3) - nf3.min_EGx );
                NE_Bound(3) = simplified_ppval( nf3.pp, nf3.min_EGx ) + abs( NE_Bound(3) );
            else if EGx(3)>nf3.max_EGx
                    NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.max_EGx ) * ( EGx(3) - nf3.max_EGx );
                    NE_Bound(3) = simplified_ppval( nf3.pp, nf3.max_EGx ) + abs( NE_Bound(3) );
                else
                    NE_Bound(3) =  simplified_ppval( nf3.pp, EGx(3) );
                end
            end

            % G5 = @(x) x*|x|/(10+|x|);
            EGx(5) = sum( y.*abs_y.*inv_p10abs_y )/T;   
            if EGx(5)<nf5.min_EGx
                NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.min_EGx ) * ( EGx(5) - nf5.min_EGx );
                NE_Bound(5) = simplified_ppval( nf5.pp, nf5.min_EGx ) + abs( NE_Bound(5) );
            else if EGx(5)>nf5.max_EGx
                    NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.max_EGx ) * ( EGx(5) - nf5.max_EGx );
                    NE_Bound(5) = simplified_ppval( nf5.pp, nf5.max_EGx ) + abs( NE_Bound(5) );
                else
                    NE_Bound(5) =  simplified_ppval( nf5.pp, EGx(5) );
                end
            end

            % G7 = @(x) x/(1+x^2);
            EGx(7) = sum( y.*inv_pabs_yy )/T;   
            if EGx(7)<nf7.min_EGx
                NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.min_EGx ) * ( EGx(7) - nf7.min_EGx );
                NE_Bound(7) = simplified_ppval( nf7.pp, nf7.min_EGx ) + abs( NE_Bound(7) );
            else if EGx(7)>nf7.max_EGx
                    NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.max_EGx ) * ( EGx(7) - nf7.max_EGx );
                    NE_Bound(7) = simplified_ppval( nf7.pp, nf7.max_EGx ) + abs( NE_Bound(7) );
                else
                    NE_Bound(7) =  simplified_ppval( nf7.pp, EGx(7) );
                end
            end

            % select the tightest bound
            [max_NE, max_i] = max( NE_Bound ); 
            negentropy1 = max_NE;
                
            y = yr2;
            % we only need to calculate these quantities once
            yy = y.*y;                              
            sign_y = sign(y);
            abs_y = sign_y.*y;
            inv_pabs_y = 1./(1+abs_y);
            inv_pabs_yy = 1./(1+yy);
            inv_p10abs_y = 1./(10+abs_y);

            % G1(x) = x^4
            EGx(1) = sum( yy.*yy )/T;  
            if EGx(1)<nf1.min_EGx
                NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.min_EGx ) * ( EGx(1) - nf1.min_EGx );
                NE_Bound(1) = simplified_ppval( nf1.pp, nf1.min_EGx ) + abs( NE_Bound(1) );
            else if EGx(1)>nf1.max_EGx
                    % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                    % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                    NE_Bound(1) = 0;
                else
                    NE_Bound(1) =  simplified_ppval( nf1.pp, EGx(1) );
                end
            end

            % G3 = @(x) abs(x)./( 1+abs(x) );
            EGx(3) = 1-sum( inv_pabs_y )/T;
            if EGx(3)<nf3.min_EGx
                NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.min_EGx ) * ( EGx(3) - nf3.min_EGx );
                NE_Bound(3) = simplified_ppval( nf3.pp, nf3.min_EGx ) + abs( NE_Bound(3) );
            else if EGx(3)>nf3.max_EGx
                    NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.max_EGx ) * ( EGx(3) - nf3.max_EGx );
                    NE_Bound(3) = simplified_ppval( nf3.pp, nf3.max_EGx ) + abs( NE_Bound(3) );
                else
                    NE_Bound(3) =  simplified_ppval( nf3.pp, EGx(3) );
                end
            end

            % G5 = @(x) x*|x|/(10+|x|);
            EGx(5) = sum( y.*abs_y.*inv_p10abs_y )/T;   
            if EGx(5)<nf5.min_EGx
                NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.min_EGx ) * ( EGx(5) - nf5.min_EGx );
                NE_Bound(5) = simplified_ppval( nf5.pp, nf5.min_EGx ) + abs( NE_Bound(5) );
            else if EGx(5)>nf5.max_EGx
                    NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.max_EGx ) * ( EGx(5) - nf5.max_EGx );
                    NE_Bound(5) = simplified_ppval( nf5.pp, nf5.max_EGx ) + abs( NE_Bound(5) );
                else
                    NE_Bound(5) =  simplified_ppval( nf5.pp, EGx(5) );
                end
            end

            % G7 = @(x) x/(1+x^2);
            EGx(7) = sum( y.*inv_pabs_yy )/T;   
            if EGx(7)<nf7.min_EGx
                NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.min_EGx ) * ( EGx(7) - nf7.min_EGx );
                NE_Bound(7) = simplified_ppval( nf7.pp, nf7.min_EGx ) + abs( NE_Bound(7) );
            else if EGx(7)>nf7.max_EGx
                    NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.max_EGx ) * ( EGx(7) - nf7.max_EGx );
                    NE_Bound(7) = simplified_ppval( nf7.pp, nf7.max_EGx ) + abs( NE_Bound(7) );
                else
                    NE_Bound(7) =  simplified_ppval( nf7.pp, EGx(7) );
                end
            end

            % select the tightest bound
            [max_NE, max_i] = max( NE_Bound ); 
            negentropy2 = max_NE;
            
            if negentropy1 + negentropy2 > max_negentropy(m) + max_negentropy(n) + eps
                if verbose
                    fprintf('\nRotating %g %g', m, n);
                end
                max_negentropy(m) = negentropy1;
                max_negentropy(n) = negentropy2;
                W(m,:) = (w1+w2)'/sqrt(2);
                W(n,:) = (w1-w2)'/sqrt(2);
                Y(m,:) = yr1;
                Y(n,:) = yr2;
                ym = yr1;
                w1 = W(m,:)';
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
Cost = zeros(max_iter_orth_refine,1);
min_cost = inf;
min_cost_queue = min_cost*ones(max_iter_orth_refine,1);
mu = 1/50;
cost_increase_counter = 0;
fastica_on = 1;
error = 0;
for iter = 1 : max_iter_orth_refine
    
    for n = 1 : N
        
        w = W(n,:)';
        y = w'*X;
        
        % evaluate the upper bound of negentropy of the nth component
        NE_Bound = zeros(K,1);
        EGx = zeros(K,1);
        
        % we only need to calculate these quantities once
        yy = y.*y;                              
        sign_y = sign(y);
        abs_y = sign_y.*y;
        inv_pabs_y = 1./(1+abs_y);
        inv_pabs_yy = 1./(1+yy);
        inv_p10abs_y = 1./(10+abs_y);
              
        % G1(x) = x^4
        EGx(1) = sum( yy.*yy )/T;  
        if EGx(1)<nf1.min_EGx
            NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.min_EGx ) * ( EGx(1) - nf1.min_EGx );
            NE_Bound(1) = simplified_ppval( nf1.pp, nf1.min_EGx ) + abs( NE_Bound(1) );
        else if EGx(1)>nf1.max_EGx
                % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                NE_Bound(1) = 0;
            else
                NE_Bound(1) =  simplified_ppval( nf1.pp, EGx(1) );
            end
        end
        
        % G3 = @(x) abs(x)./( 1+abs(x) );
        EGx(3) = 1-sum( inv_pabs_y )/T;
        if EGx(3)<nf3.min_EGx
            NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.min_EGx ) * ( EGx(3) - nf3.min_EGx );
            NE_Bound(3) = simplified_ppval( nf3.pp, nf3.min_EGx ) + abs( NE_Bound(3) );
        else if EGx(3)>nf3.max_EGx
                NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.max_EGx ) * ( EGx(3) - nf3.max_EGx );
                NE_Bound(3) = simplified_ppval( nf3.pp, nf3.max_EGx ) + abs( NE_Bound(3) );
            else
                NE_Bound(3) =  simplified_ppval( nf3.pp, EGx(3) );
            end
        end
        
        % G5 = @(x) x*|x|/(10+|x|);
        EGx(5) = sum( y.*abs_y.*inv_p10abs_y )/T;   
        if EGx(5)<nf5.min_EGx
            NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.min_EGx ) * ( EGx(5) - nf5.min_EGx );
            NE_Bound(5) = simplified_ppval( nf5.pp, nf5.min_EGx ) + abs( NE_Bound(5) );
        else if EGx(5)>nf5.max_EGx
                NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.max_EGx ) * ( EGx(5) - nf5.max_EGx );
                NE_Bound(5) = simplified_ppval( nf5.pp, nf5.max_EGx ) + abs( NE_Bound(5) );
            else
                NE_Bound(5) =  simplified_ppval( nf5.pp, EGx(5) );
            end
        end
        
        % G7 = @(x) x/(1+x^2);
        EGx(7) = sum( y.*inv_pabs_yy )/T;   
        if EGx(7)<nf7.min_EGx
            NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.min_EGx ) * ( EGx(7) - nf7.min_EGx );
            NE_Bound(7) = simplified_ppval( nf7.pp, nf7.min_EGx ) + abs( NE_Bound(7) );
        else if EGx(7)>nf7.max_EGx
                NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.max_EGx ) * ( EGx(7) - nf7.max_EGx );
                NE_Bound(7) = simplified_ppval( nf7.pp, nf7.max_EGx ) + abs( NE_Bound(7) );
            else
                NE_Bound(7) =  simplified_ppval( nf7.pp, EGx(7) );
            end
        end
               
        % select the tightest bound
        [max_NE, max_i] = max( NE_Bound );  
        negentropy_array(n) = max_NE;
        Cost(iter) = Cost(iter) - max_NE;
            
            % Perform orthogonal ICA
            switch max_i
                
                case 1
                    
                    % G1(x) = x^4
                    % g1(x) = 4x^3
                    % g1' =  12x^2
                    grad = X*( 4*y.*yy )'/T;
                    if fastica_on 
                        Edgx = 12; 
                    else
                        vEGx = 2*( EGx(1) > nf1.critical_point ) - 1;
                    end
                    
                case 3
                    
                    % G3 = @(x) abs(x)./( 1+abs(x) );
                    % g3 = @(x) sign(x)./( 1+abs(x) ).^2;
                    % g3' = -2/(1+abs(x))^3
                    grad = X*( sign_y.*inv_pabs_y.*inv_pabs_y )'/T;      %grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
                    if fastica_on 
                        Edgx = sum(-2*inv_pabs_y.*inv_pabs_y.*inv_pabs_y)/T; % Edgx = sum(-2*inv_pabs_y.^3)/T; 
                    else
                        vEGx = 2*( EGx(3) > nf3.critical_point ) - 1;
                    end
                    
                case 5
                    
                    % G5 = @(x) x*|x|/(10+|x|);
                    % g5 = @(x) |x|*(20+|x|)/(10+|x|)^2;      
                    % g5 = 200*sign(x) / (10+|x|)^3
                    grad = X*( abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'/T;
                    if fastica_on 
                        Edgx = sum( 200*sign_y.*inv_p10abs_y.*inv_p10abs_y.*inv_p10abs_y )/T;    %Edgx = sum( 200*sign_y.*inv_p10abs_y.^3 )/T; 
                    else
                        vEGx = 2*( EGx(5) > nf5.critical_point ) - 1;
                    end
                    
                case 7
                    
                    % G7 = @(x) x/(1+x^2)
                    % g7 = @(x) (1-x^2)/(1+x^2)^2
                    % g7' = 2x(x^2-3)/(1+x^2)^3
                    grad = X*( (1-yy).*inv_pabs_yy.^2 )'/T;
                    if fastica_on 
                        Edgx = sum( 2*y.*(yy-3).*inv_pabs_yy.*inv_pabs_yy.*inv_pabs_yy )/T;    %Edgx = sum( 2*y.*(yy-3).*inv_pabs_yy.^3 )/T; 
                    else
                        vEGx = 2*( EGx(7) > nf7.critical_point ) - 1;
                    end
                                        
                otherwise
                    ;
            end
            
            if fastica_on
                w1 = grad - Edgx*w;
            else  
                grad = vEGx*grad;
                grad = grad - (w'*grad)*w;
                grad = grad / norm(grad);
                w1 = w + mu*grad;
            end
              
            W(n,:) = w1';
 
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
        if cost_increase_counter>=max_cost_increase_number || 1-min(abs(diag(W*last_W'))) < tolerance
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
    subplot(1,4,3);
    plot([1:iter_fastica],Cost([1:iter_fastica]), 'r');
    hold on; plot([iter_fastica+1:iter],Cost([iter_fastica+1:iter]), 'b');
    xlabel('Number of iterations')
    ylabel('Cost')
    title('(b) Refine orthogonal ICA')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
if show_cost
    subplot(1,4,3);
    title('(b) No saddles, no Refinement')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sort all the component
[max_negentropy, index_sort] = sort(max_negentropy, 'descend');
W = W(index_sort, :);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%           Part IV: non-orth ICA
%   fixed small step size for refinement, gradient search
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if verbose
    fprintf('\nNonorthogonal ICA stage.');
end
last_W = W;
best_W = W;
Cost = zeros(max_iter_nonorth,1);
% min_cost = inf;
min_cost_queue = min_cost*ones(max_iter_nonorth,1);
error = inf;
mu = 1/25;    %step size
min_mu = 1/200;
max_cost_increase_number = 3;
cost_increase_counter = 0;
for iter = 1 : max_iter_nonorth
    
    Cost(iter) = -log(abs(det(W)));
    % CONSTRAINT
   
    %W = W(randperm(N),:);
    for n = 1 : N
%          if n <= num_guess
%             r_n_c = guess_mat(:,n);
%          end
        % if N is too large, we adopt a smarter way to calculate h_n to
        % reduce the computational load
        if N>7
            
            if n==1
                Wn = W(2:N,:);
                inv_Q = inv( Wn*Wn' );
            else
                n_last = n-1;
                Wn_last = [ W(1:n_last-1,:);W(n_last+1:N,:) ];
                w_current = W(n,:)';
                w_last = W(n_last,:)';
                c = Wn_last*( w_last - w_current );
                c(n_last) = 0.5*( w_last'*w_last - w_current'*w_current );
                e_last = zeros(N-1,1);
                e_last(n_last) = 1;
            
                temp1 = inv_Q*c;
                temp2 = inv_Q(:,n_last);
                inv_Q_plus = inv_Q - temp1*temp2'/(1+temp1(n_last));
            
                temp1 = inv_Q_plus'*c;
                temp2 = inv_Q_plus(:,n_last);
                inv_Q = inv_Q_plus - temp2*temp1'/(1+c'*temp2);
                % inv_Q is Hermitian
                inv_Q = (inv_Q+inv_Q')/2;
            end

            temp1 = randn(N, 1);
            W_n = [ W(1:n-1,:);W(n+1:N,:) ];
            h = temp1 - W_n'*inv_Q*W_n*temp1;
            
        else
            temp1 = randn(N, 1);
            temp2 = [W(1:n-1,:); W(n+1:N,:)];
            h = temp1 - temp2'*inv(temp2*temp2')*temp2*temp1;
        end
        
        w = W(n,:)';
        y = w'*X;
        
        % evaluate the upper bound of negentropy of the n-th component
        NE_Bound = zeros(K,1);
        EGx = zeros(K,1);
        
        % we only need to calculate these quantities once
        yy = y.*y;                              
        sign_y = sign(y);
        abs_y = sign_y.*y;
        inv_pabs_y = 1./(1+abs_y);
        inv_pabs_yy = 1./(1+yy);
        inv_p10abs_y = 1./(10+abs_y);
        
        % G1(x) = x^4
        EGx(1) = sum( yy.*yy )/T;  
        if EGx(1)<nf1.min_EGx
            NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.min_EGx ) * ( EGx(1) - nf1.min_EGx );
            NE_Bound(1) = simplified_ppval( nf1.pp, nf1.min_EGx ) + abs( NE_Bound(1) );
        else if EGx(1)>nf1.max_EGx
                % NE_Bound(1) =  simplified_ppval( nf1.pp_slope, nf1.max_EGx ) * ( EGx(1) - nf1.max_EGx );
                % NE_Bound(1) = simplified_ppval( nf1.pp, nf1.max_EGx ) + abs( NE_Bound(1) );
                NE_Bound(1) = 0;
            else
                NE_Bound(1) =  simplified_ppval( nf1.pp, EGx(1) );
            end
        end
        
        % G3 = @(x) abs(x)./( 1+abs(x) );
        EGx(3) = 1-sum( inv_pabs_y )/T;
        if EGx(3)<nf3.min_EGx
            NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.min_EGx ) * ( EGx(3) - nf3.min_EGx );
            NE_Bound(3) = simplified_ppval( nf3.pp, nf3.min_EGx ) + abs( NE_Bound(3) );
        else if EGx(3)>nf3.max_EGx
                NE_Bound(3) =  simplified_ppval( nf3.pp_slope, nf3.max_EGx ) * ( EGx(3) - nf3.max_EGx );
                NE_Bound(3) = simplified_ppval( nf3.pp, nf3.max_EGx ) + abs( NE_Bound(3) );
            else
                NE_Bound(3) =  simplified_ppval( nf3.pp, EGx(3) );
            end
        end
        
        % G5 = @(x) x*|x|/(10+|x|);
        EGx(5) = sum( y.*abs_y.*inv_p10abs_y )/T;   
        if EGx(5)<nf5.min_EGx
            NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.min_EGx ) * ( EGx(5) - nf5.min_EGx );
            NE_Bound(5) = simplified_ppval( nf5.pp, nf5.min_EGx ) + abs( NE_Bound(5) );
        else if EGx(5)>nf5.max_EGx
                NE_Bound(5) =  simplified_ppval( nf5.pp_slope, nf5.max_EGx ) * ( EGx(5) - nf5.max_EGx );
                NE_Bound(5) = simplified_ppval( nf5.pp, nf5.max_EGx ) + abs( NE_Bound(5) );
            else
                NE_Bound(5) =  simplified_ppval( nf5.pp, EGx(5) );
            end
        end
        
        % G7 = @(x) x/(1+x^2);
        EGx(7) = sum( y.*inv_pabs_yy )/T;   
        if EGx(7)<nf7.min_EGx
            NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.min_EGx ) * ( EGx(7) - nf7.min_EGx );
            NE_Bound(7) = simplified_ppval( nf7.pp, nf7.min_EGx ) + abs( NE_Bound(7) );
        else if EGx(7)>nf7.max_EGx
                NE_Bound(7) =  simplified_ppval( nf7.pp_slope, nf7.max_EGx ) * ( EGx(7) - nf7.max_EGx );
                NE_Bound(7) = simplified_ppval( nf7.pp, nf7.max_EGx ) + abs( NE_Bound(7) );
            else
                NE_Bound(7) =  simplified_ppval( nf7.pp, EGx(7) );
            end
        end
                       
        [max_NE, max_i] = max( NE_Bound );      % select the tightest bound
%         Cost(iter) = Cost(iter) - max_NE;
  
         % CONSTRAINT         
         if n <= num_guess
%              if constraint_type == 0
%                  e_pair = corrcoef(V_1(:,:,k)*W(n,:,k)',r_n_c); % correlation of w and r_n,W(n,:,k)'\V(:,:,k)
%                  dis_wr=abs(e_pair(1,2));
%              else
                 r_n_c = guess_mat(:,n);
                 e_pair = corrcoef(y',r_n_c); %correlation of y and r_n
                 dis_wr=abs(e_pair(1,2));
%              end
             mu_old(n) = mu_c;
             mu_new(n)=max(0,mu_c + gam*(rho-dis_wr));
             mu_c=sign(e_pair(1,2))*mu_new(n);

             %Modify r_i: whitening
%              if constraint_type == 0 % constrain w
%                  r_n_mu=VV(:,:,k)*r_n_c; % If r_n is for w, it should be whiten,r_n\(V(:,:,k)')
%                  r_n_mu=r_n_mu/norm(r_n_mu);
%                  w1 = w + mu*grad - mu_c(n).*r_n_mu; % A column vector
%              else
                 r_n_c=r_n_c/norm(r_n_c);
%                  w1 = w + mu*grad - mu_c.*(X*r_n_c); % A column vector
%              end 
         end
%         Cost(iter) = Cost(iter) - max_NE - (sum(sum(mu_new(n).^2 - mu_old(n).^2)))/(2*gam);

                    
            switch max_i
                
                case 1
                    
%                     % G1(x) = x^4
%                     % g1(x) = 4x^3
%                     % g1' =  12x^2
%                     vEGx = 2*( EGx(1) > nf1.critical_point ) - 1;
%                     grad = X*( 4*y.*yy )'/T;  
                    EGx(1) = max( min( EGx(1), nf1.max_EGx ), nf1.min_EGx );
                    grad = h/(h'*w) + X*( 4*y.*yy )'*simplified_ppval( nf1.pp_slope, EGx(1) ) /T;      % gradient
                    
                case 3
                    
%                     % G3 = @(x) abs(x)./( 1+abs(x) );
%                     % g3 = @(x) sign(x)./( 1+abs(x) ).^2;
%                     % g3' = -2/(1+abs(x))^3
%                     vEGx = 2*( EGx(3) > nf3.critical_point ) - 1;
%                     grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
                    EGx(3) = max( min( EGx(3), nf3.max_EGx ), nf3.min_EGx );
                    grad = h/(h'*w) + X *( sign_y.*inv_pabs_y.^2 )'*simplified_ppval( nf3.pp_slope, EGx(3) )/T;
                    
                case 5
                    
%                     % G5 = @(x) x*|x|/(10+|x|);
%                     % g5 = @(x) |x|*(20+|x|)/(10+|x|)^2;      
%                     % g5 = 200*sign(x) / (10+|x|)^3
%                     vEGx = 2*( EGx(5) > nf5.critical_point ) - 1;
%                     grad = X*( abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'/T;
                    EGx(5) = max( min( EGx(5), nf5.max_EGx ), nf5.min_EGx );
                    grad = h/(h'*w) + X*( abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'*simplified_ppval( nf5.pp_slope, EGx(5) ) /T;
                                        
                case 7
                    
%                     % G7 = @(x) x/(1+x^2)
%                     % g7 = @(x) (1-x^2)/(1+x^2)^2
%                     vEGx = 2*( EGx(7) > nf7.critical_point ) - 1;
%                     grad = X*( (1-yy).*inv_pabs_yy.^2 )'/T; 
                    EGx(7) = max( min( EGx(7), nf7.max_EGx ), nf7.min_EGx );
                    grad = h/(h'*w) + X*( (1-yy).*inv_pabs_yy.^2 )'*simplified_ppval( nf7.pp_slope, EGx(7) ) /T;   %this is the negative value of eq14
                                                           
                otherwise
                    ;
            end
         grad = grad - (w'*grad)*w;
         grad = grad / norm(grad);
       
         
           w1 = w + mu*grad - mu_c.*(X*r_n_c); % A column vector

            w1 = w1 / norm(w1);
            W(n,:) = w1';
            Cost(iter) = Cost(iter) - max_NE - (sum(mu_new(n).^2 - mu_old(n).^2))/(2*gam);

    end
%     Cost(iter) = Cost(iter) - max_NE;

%     Cost(iter) = Cost(iter) - (sum(sum(mu_new.^2 - mu_old.^2)))/(2*gam);

    
    if Cost(iter) < min_cost
        cost_increase_counter = 0;
        min_cost = Cost(iter);
        best_W = last_W;
    else
        cost_increase_counter = cost_increase_counter + 1;
    end
    min_cost_queue(iter) = min_cost;
    
    if cost_increase_counter > max_cost_increase_number
        if mu > min_mu
            cost_increase_counter = 0;
            W = best_W;
            last_W = W;
            mu = mu/2;
            continue;
        else
            break;
        end
    else
        last_W = W;
    end

end
W = best_W;
W = W*P;

if show_cost
    subplot(1,4,4);
    plot([1:iter],Cost([1:iter]));
    xlabel('Number of iterations')
    ylabel('Cost')
    title('(c) Nonorthogonal ICA')
end
% [EOF ICA]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this subfunction is faster than ppval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function v = simplified_ppval(pp,xs)
% a simplified version of ppval

b = pp.breaks;
c = pp.coefs;
l = pp.pieces;
k = 4;  % k = pp.order;
dd = 1; % dd = pp.dim;

% find index
index=0;
middle_index=0;
if xs>b(l)
    index = l;
else if xs<b(2)
        index = 1;
    else 
        
        low_index = 1;
        high_index = l;
        while 1
            
            middle_index = round( 0.6*low_index + 0.4*high_index );
            if b( middle_index ) > xs
                high_index = middle_index;
            else
                low_index = middle_index;
            end
            
            if low_index == high_index-1
                index = low_index;
                break;
            end
            
        end
        
    end
end
    

% now go to local coordinates ...
xs = xs-b(index);

% ... and apply nested multiplication:
   v = c(index,1);
   for i=2:k
      v = xs*v + c(index,i);
   end
   
   
function W=symdecor(M)
%fast symmetric orthogonalization  
[V D]=eig(M*M');
W=(V.*(ones(size(M,2),1)*(1./sqrt(diag(D)'))))*V'*M;


function [X, P] = pre_processing( X )
% pre-processing program
[N,T] = size(X);
% remove DC
Xmean=mean(X,2);
X = X - Xmean*ones(1,T);    

% spatio pre-whitening 1 
R = X*X'/T;                 
P = inv_sqrtmH(R);  %P = inv(sqrtm(R));
X = P*X;



function A = inv_sqrtmH( B )
% 
[V,D] = eig(B);
d = diag(D);
d = 1./sqrt(d);
A = V*diag(d)*V';

function properties = getopt(properties,varargin)
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
for ii=1:length(varargin)
   arg = varargin{ii};
   if isempty(TargetField)
      if ~ischar(arg)
         error('Property names must be character strings');
      end
      %f = find(strcmp(prop_names, arg));
      if isempty(find(strcmp(prop_names, arg),1)) %length(f) == 0
         error('%s ',['invalid property ''',arg,'''; must be one of:'],prop_names{:});
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
