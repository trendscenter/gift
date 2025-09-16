function W = icatb_fbss(X, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BSS by fully exploiting the non-Gausianity and correlation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs:
% X:        mixtures;
% p:        filter length, p is a positive integer
% Output:
% W:        demixing matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference:
%
% Xi-Lin Li, Tulay Adali, "Blind spatiotemporal separation of second and/or
% higher-order correlated sources by entropy rate minimization,"
% IEEE International Conference on Acoustics, Speech and Signal Processing 2010.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Coded by Xilin Li, lixilin@umbc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     Part 0: pre-processing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global nf1 nf2 nf3 nf4 nf5 nf6 nf7 nf8
show_cost = 0;

[N,T] = size(X);
[X, P] = pre_processing( X );

p = [];
for nV = 1:length(varargin)
    if (strcmpi(varargin{nV}, 'filter_length'))
        p = varargin{nV + 1};
    end
end

if (isempty(p))
    p = min(11, T/50);
end

% % read the input arguments and initialization
% if nargin == 1
%     p = min(11, T/50);
% end

W = icatb_ica_ebm(X);     % use ICA-EBM to provide the initial guess
if p == 1               % only exploiting the non-Gaussianity for separation
    W = W*P;
    return;
end

% Load 8 measuring functions. But we only use 4 of them
K = 8;
load icatb_nf_table;

% prediction coeffs
a = zeros(p, N);
for n = 1 : N
    a(floor((p+1)/2), n) = 1;
end

Rz = zeros(N,N,N);
temp5 = zeros(T,N);
Z = zeros(N,T,N);

calculate_cos_sin_mtx( p ); % prepare the data used in integral

last_W = W;
best_W = W;
best_a = a;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Part I: Here begins our work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for stochastic_search = 1 : -1: 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if stochastic_search
        mu = 1/5;               % use stochastic search to avoid local convergence
        max_cost_increase = 5;
        max_iter_north = 500;
        tolerance = 1e-3;
    else
        mu = 1/50;
        max_cost_increase = 3;
        max_iter_north = 200;
        tolerance = 1e-5;
    end
    cost_increase_counter = 0;
    W = best_W;
    a = best_a;
    last_W = best_W;
    Cost = zeros(max_iter_north,1);
    min_cost = inf;
    min_cost_queue = min_cost*ones(max_iter_north,1);
    negentropy_array = zeros(N,1);
    
    for iter = 1 : max_iter_north
        
        if stochastic_search
            % estimate the AR coeffs
            Y = W*X;
            for n = 1 : N
                
                %             if min(abs(1-abs(roots(a(:,n))))) < 2e-2
                %                 continue;
                %             end
                
                if mod(iter,6)==1 || iter<=5
                    
                    [a1, min_ere1] = lfc(Y(n,:), p, 'unknown', []);         % try a new one
                    [a2, min_ere2] = lfc(Y(n,:), p, [], a(:, n));           % refine
                    
                    % choose the best model
                    min_ere = inf;
                    if min_ere > min_ere1
                        min_ere = min_ere1;
                        a(:, n) = a1;
                    end
                    if min_ere > min_ere2
                        min_ere = min_ere2;
                        a(:, n) = a2;
                    end
                    
                else if mod(iter,6)==4
                        
                        a(:, n) = lfc(Y(n,:), p, [], a(:, n));  % refine. we cannot afford to do it every iteration
                        
                    end
                end
                
                temp5 = filter(a(:,n), 1, X');
                Rz(:,:,n) = temp5'*temp5/T;
                Z(:,:,n) = temp5';
            end
        end
        
        Cost(iter) = -log(abs(det(W)));
        
        % estimate W
        for n = 1 : N
            
            temp1 = randn(N, 1);
            temp2 = [W(1:n-1,:); W(n+1:N,:)];
            h = temp1 - temp2'*((temp2*temp2')\temp2)*temp1;            
            v = W(n,:)';
            sigma2 = v'*Rz(:,:,n)*v;
            Cost(iter) = Cost(iter) + log(sigma2)/2;
            v = v/sqrt(sigma2);
            y = v'*Z(:,:,n);           % this y is the prediction error
            
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
            negentropy_array(n) = max_NE;
            Cost(iter) = Cost(iter) - max_NE;
            
            if stochastic_search
                weight = rand(1,T);
            else
                weight = ones(1,T);
            end
            
            switch max_i
                
                case 1
                    
                    % G1(x) = x^4
                    % g1(x) = 4*x^3
                    EGx(1) = max( min( EGx(1), nf1.max_EGx ), nf1.min_EGx );
                    grad = h/(h'*v) + Z(:,:,n)*( 4*weight.*y.*yy )'*simplified_ppval( nf1.pp_slope, EGx(1) ) /sum(weight);     % gradient
                    
                case 3
                    
                    % G3 = @(x) abs(x)./( 1+abs(x) );
                    % g3 = @(x) sign(x)./( 1+abs(x) ).^2;
                    EGx(3) = max( min( EGx(3), nf3.max_EGx ), nf3.min_EGx );
                    grad = h/(h'*v) + Z(:,:,n) *( weight.*sign_y.*inv_pabs_y.^2 )'*simplified_ppval( nf3.pp_slope, EGx(3) )/sum(weight);
                    
                case 5
                    
                    % G5 = x*|x|/(10+|x|);
                    % g5 = |x|*(20+|x|)/(10+|x|)^2;
                    EGx(5) = max( min( EGx(5), nf5.max_EGx ), nf5.min_EGx );
                    grad = h/(h'*v) + Z(:,:,n)*( weight.*abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'*simplified_ppval( nf5.pp_slope, EGx(5) ) /sum(weight);
                    
                case 7
                    
                    % G7 = @(x) x/(1+x^2)
                    % g7 = @(x) (1-x^2)/(1+x^2)^2
                    EGx(7) = max( min( EGx(7), nf7.max_EGx ), nf7.min_EGx );
                    grad = h/(h'*v) + Z(:,:,n)*( weight.*(1-yy).*inv_pabs_yy.^2 )'*simplified_ppval( nf7.pp_slope, EGx(7) ) /sum(weight);
                    
                otherwise
                    ;
            end
            
            cnstd = Rz(:,:,n)*v;    % constraint direction
            grad = grad - cnstd'*grad*cnstd/(cnstd'*cnstd); % projected gradient
            grad = inv_sqrtmH( Rz(:,:,n) )*grad;    %grad = inv(sqrtm( Rz(:,:,n)+max(diag(Rz(:,:,n)))*eye(N)/100 ))*grad;  % this whitening operation is to speed up convergence
            grad = grad / sqrt(grad'*Rz(:,:,n)*grad);   % normalized gradient
            v = v + mu*grad;
            v = v/sqrt(v'*Rz(:,:,n)*v);
            W(n,:) = v';
            
        end
        
        if Cost(iter) < min_cost
            cost_increase_counter = 0;
            min_cost = Cost(iter);
            best_W = last_W;
            best_a = a;
            max_negentropy = negentropy_array;
        else
            cost_increase_counter = cost_increase_counter + 1;
        end
        min_cost_queue(iter) = min_cost;
        
        if cost_increase_counter > max_cost_increase
            if stochastic_search
                W1 = W;
                last_W1 = last_W;
                for n = 1 : N
                    W1(n,:) = W1(n,:)/norm(W1(n,:));
                    last_W1(n,:) = last_W1(n,:)/norm(last_W1(n,:));
                end
                if 1-min(abs(diag(W1*last_W1'))) < tolerance
                    break;
                else
                    mu = mu/2;
                    W = best_W;
                    last_W = best_W;
                    a = best_a;
                    cost_increase_counter = 0;
                    continue;
                end
            else
                W1 = W;
                last_W1 = last_W;
                for n = 1 : N
                    W1(n,:) = W1(n,:)/norm(W1(n,:));
                    last_W1(n,:) = last_W1(n,:)/norm(last_W1(n,:));
                end
                if 1-min(abs(diag(W1*last_W1'))) < tolerance
                    break;
                else
                    mu = mu/2;
                    W = best_W;
                    last_W = best_W;
                    a = best_a;
                    cost_increase_counter = 0;
                    continue;
                end
            end
        end
        
        last_W = W;
        
    end
    W = best_W;
    
    if show_cost
        if stochastic_search
            figure; subplot(1,2,1);
        else
            subplot(1,2,2);
        end
        plot([1:iter],Cost([1:iter]));
        xlabel('Number of iterations')
        ylabel('Cost')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end     % end of work
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % sort all component
% for n = 1 : N
%     max_negentropy(n) = max_negentropy(n) + log(norm(W(n,:)));
% end
% [max_negentropy, index_sort] = sort(max_negentropy, 'descend');
% W = W(index_sort, :);

W = W*P;

% [EOF FBSS]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [X, P] = pre_processing( X )
% pre-processing program
[N,T] = size(X);
% remove DC
Xmean=mean(X,2);
X = X - Xmean*ones(1,T);

% spatio pre-whitening 1
R = X*X'/T;
P1 = inv_sqrtmH(R);  %P = inv(sqrtm(R));
X = P1*X;
% temporal pre-filtering for colored signals
q = 3;
r = zeros(q,1);
for p = 0 : q-1
    r(p+1) = trace(X(:,1:T-p)  *X(:,p+1:T)')/T/N;
end
%af = inv(toeplitz(r(1:q-1)))*conj(r(2:q));
af = toeplitz(r(1:q-1)) \ conj(r(2:q));
for n = 1 : N
    X(n,:) = filter( [1;-af],1,X(n,:) );
end
% spatio pre-whitening 2
R = X*X'/T;
P2 = inv_sqrtmH(R);  %P = inv(sqrtm(R));
X = P2*X;
P = P2*P1;



function A = inv_sqrtmH( B )
%
[V,D] = eig(B);
d = diag(D);
d = 1./sqrt(d);
A = V*diag(d)*V';


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
            
            middle_index = round( (low_index + high_index)/2 );
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
    v = xs.*v + c(index,i);
end










function calculate_cos_sin_mtx( p )
% prepare the cos and sin matrix for integral calculation

global cosmtx
global sinmtx
global Simpson_c

% sampling omega from 0 to pi
n = 10*p;           % 10*length(a) samples are enough
h = pi/n;
w = ([0:n])*h;

cosmtx = zeros(p, n+1);
sinmtx = zeros(p, n+1);
for q = 1 : p
    cosmtx(q, :) = cos((q-1)*w);
    sinmtx(q, :) = sin((q-1)*w);
end

% c is the vector like [1,4,2,4,2,...,2,4,1] used in Simpson's Rule
Simpson_c = zeros(n+1,1);
Simpson_c(1:2:n+1) = 2;
Simpson_c(2:2:n) = 4;
Simpson_c(1) = 1;
Simpson_c(n+1) = 1;











function [b, G] = cnstd_and_gain (a)
% return constraint direction used for calculating projected gradient and Gain of filter a

global cosmtx
global sinmtx
global Simpson_c

p = length(a);

% calculate the integral

% sampling omega from 0 to pi
n = 10*p;           % 10*length(a) samples are enough
h = pi/n;
w = ([0:n])*h;

% calculate |A(w)|^2
Awr = zeros(1,n+1);        % Re[A(w)]
Awi = zeros(1,n+1);        % Im[A(w)]
for q = 1 : p
    Awr = Awr + a(q)*cosmtx(q, :);
    Awi = Awi + a(q)*sinmtx(q, :);
end
Aw2 = 10*eps + Awr.^2 + Awi.^2;       % add 100*eps to prevent log of zero

% calculate the vector
v = zeros(p+1,n+1);
inv_Aw2 = 1./Aw2;
for q = 1 : p
    %v(q,:) = cos((q-1)*w)./Aw2;
    v(q,:) = cosmtx(q, :).*inv_Aw2;
end
v(p+1,:) = log(Aw2)/pi;

% this is the integral
u = h*v*Simpson_c/3;

b = toeplitz(u(1:p))*a; % this is the vector
G = u(p+1);             % this is the gain







function [a, min_cost] = lfc(x, p, choice, a0)
% Return the linear filtering coefficient (LFC) with length p for entropy
% rate estimation, and the estimated entropy rate
% Inputs
% p is the filter length
% 'choice' can be 'sub', 'super', or 'unknown'
% a0 is the intial guess
% Outputs
% a is the filter coefficients
% min_cost is the entropy rate estimation
global nf1 nf2 nf3 nf4 nf5 nf6 nf7 nf8
tolerance = 1e-4;
T = length(x);

X0 = icatb_convmtx(x,p);
X = X0(:, 1:T);         % if use the whole convmtx, outliers arise, so remove the tail
% remove DC
Xmean=mean(X,2);
X = X - Xmean*ones(1,T);
% pre-whitening
R = X*X'/T;
[V,D] = eig(R);
d = diag(D);
d(d<10*eps) = 10*eps;
P = V*diag(1./sqrt( d ))*V';
X = P*X;

if isempty(a0)
    % use SEA to provide the initial guess
    
    switch choice
        
        case 'sub'
            a = minkurt( X );
            
        case 'super'
            a = maxkurt( X );
            
        otherwise
            % a = zeros(p,1);
            % a(floor((p+1)/2)) = 1;
            a = randn(p,1);     % since it will run multiple times in my ICA code, I prefer a random initialization
            a = a / norm(a);
            last_a = a;
            for iter = 1 : 100
                y = a'*X;
                a = X*(y.^3)'/T - 3*a;
                a = a / norm(a);
                if 1-abs(a'*last_a) < tolerance
                    break;
                else
                    last_a = a;
                end
            end
            
    end
    
else
    %a = inv(P)*a0;
    a = P \ a0;
end

% Load 8 measuring functions.
K = 8;
% load nf_table;
% load nf2_bounded;   % unbounded function |x| is replaced with min(|x|,10), a bounded one

min_cost = inf;
best_a = a;
last_a = a;
min_mu = 1/128;
if isempty(a0)
    max_iter = 100;
    mu = 4*min_mu;
else
    max_iter = 100;
    mu = 16*min_mu;
end
cost_increase_counter = 0;
Cost = zeros(max_iter,1);
for iter = 1 : max_iter
    
    a_original = P*a;       % a_original is the a before whitening
    [b_orginal, G_orginal] = cnstd_and_gain (a_original);
    a = a*exp(-G_orginal/2);
    b = P*b_orginal;
    
    y = a'*X;
    sigma2 = a'*a;
    y = y/sqrt(sigma2);     % y is always the normalized y here
    Cost(iter) = 0.5*log(2*pi*sigma2)+0.5;
    
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
    Cost(iter) = Cost(iter) - max_NE;
    
    last_a = a;
    
    if Cost(iter) < min_cost
        cost_increase_counter = 0;
        min_cost = Cost(iter);
        best_a = a;
    else
        cost_increase_counter = cost_increase_counter + 1;
    end
    
    if cost_increase_counter > 0
        if mu > min_mu
            mu = mu/2;
            cost_increase_counter = 0;
            a = best_a;
            last_a = best_a;
            continue;
        else
            break;
        end
    end
    
    grad = a/sigma2;
    
    switch max_i
        
        case 1
            
            %                     % G1(x) = x^4
            %                     % g1(x) = 4x^3
            %                     % g1' =  12x^2
            %                     vEGx = 2*( EGx(1) > nf1.critical_point ) - 1;
            %                     grad = X*( 4*y.*yy )'/T;
            EGx(1) = max( min( EGx(1), nf1.max_EGx ), nf1.min_EGx );
            grad = grad - X*( 4*y.*yy )'*simplified_ppval( nf1.pp_slope, EGx(1) ) /T/sqrt(sigma2);
            grad = grad + sum(4*y.*yy.*y)*simplified_ppval( nf1.pp_slope, EGx(1) )*a/T/sigma2;
            
        case 3
            
            %                     % G3 = @(x) abs(x)./( 1+abs(x) );
            %                     % g3 = @(x) sign(x)./( 1+abs(x) ).^2;
            %                     % g3' = -2/(1+abs(x))^3
            %                     vEGx = 2*( EGx(3) > nf3.critical_point ) - 1;
            %                     grad = X*( sign_y.*inv_pabs_y.^2 )'/T;
            EGx(3) = max( min( EGx(3), nf3.max_EGx ), nf3.min_EGx );
            grad = grad - X *( sign_y.*inv_pabs_y.^2 )'*simplified_ppval( nf3.pp_slope, EGx(3) )/T/sqrt(sigma2);
            grad = grad + sum( sign_y.*inv_pabs_y.^2.*y )*simplified_ppval( nf3.pp_slope, EGx(3) )*a/T/sigma2;
            
        case 5
            
            %                     % G5 = @(x) x*|x|/(10+|x|);
            %                     % g5 = @(x) |x|*(20+|x|)/(10+|x|)^2;
            %                     % g5 = 200*sign(x) / (10+|x|)^3
            %                     vEGx = 2*( EGx(5) > nf5.critical_point ) - 1;
            %                     grad = X*( abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'/T;
            EGx(5) = max( min( EGx(5), nf5.max_EGx ), nf5.min_EGx );
            grad = grad - X*( abs_y.*(20+abs_y).*inv_p10abs_y.^2 )'*simplified_ppval( nf5.pp_slope, EGx(5) ) /T/sqrt(sigma2);
            grad = grad + sum( abs_y.*(20+abs_y).*inv_p10abs_y.^2.*y )*simplified_ppval( nf5.pp_slope, EGx(5) )*a /T/sigma2;
            
        case 7
            
            %                     % G7 = @(x) x/(1+x^2)
            %                     % g7 = @(x) (1-x^2)/(1+x^2)^2
            %                     vEGx = 2*( EGx(7) > nf7.critical_point ) - 1;
            %                     grad = X*( (1-yy).*inv_pabs_yy.^2 )'/T;
            EGx(7) = max( min( EGx(7), nf7.max_EGx ), nf7.min_EGx );
            grad = grad - X*( (1-yy).*inv_pabs_yy.^2 )'*simplified_ppval( nf7.pp_slope, EGx(7) ) /T/sqrt(sigma2);
            grad = grad + sum( (1-yy).*inv_pabs_yy.^2.*y )*simplified_ppval( nf7.pp_slope, EGx(7) )*a /T/sigma2;
            
        otherwise
            ;
    end
    
    grad = grad - grad'*b*b/(b'*b);
    grad = sqrt(sigma2) * grad / norm(grad);
    a = a - mu*grad;
    
end
a = best_a;
a = P*a;

% figure;
% plot([1:iter], Cost(1:iter), '-o')
% pause(1/2)
%
% figure
% hist(best_a'*X,100)