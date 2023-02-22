function [comp_est, mdl, aic, kic] = icatb_estimate_dimension(files, maskvec, precisionType, dim_est_opts)
% function [comp_est, mdl, aic, kic] = icatb_estimate_dimension(files, maskvec)
% Select the order of the multivariate data using Information theoretic
% criteria with the option of sample dependence correction;
% Inputs:
%
% 1. files - fullfile path of data files as a character array or
% numeric array of size xdim by ydim by zdim by tdim.
% 2. maskvec - Mask file in analyze or Nifti format. Mask must have the
% same dimensions as the data or mask vector must be in logical vector format of
% of length xdim*ydim*zdim.
% 3. precisionType - Load data as double or single
% 4. dim_est_opts - Dimensionality options. Options to select method and if
% you specified method = 2 set fwhm variable as well.
%
% Output:
% 1. comp_est:  estimated order using MDL
% 2. mdl - MDL vector
% 3. aic - AIC vector
% 4. kic - KIC vector

% Please cite the following paper if you use this code for publication
% Y.-O. Li, T. Adali and V. D. Calhoun, "Estimating the number of independent
% components for fMRI data," Human Brain Mapping, vol. 28, pp. 1251--66, 2007


%COPYRIGHT NOTICE
%This function is a part of GIFT software library
%Copyright (C) 2003-2009
%
%This program is free software; you can redistribute it and/or
%modify it under the terms of the GNU General Public License
%as published by the Free Software Foundation; either version 2
%of the License, or any later version.
%
%This program is distributed in the hope that it will be useful,
%but WITHOUT ANY WARRANTY; without even the implied warranty of
%MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%GNU General Public License for more details.
%
%You should have received a copy of the GNU General Public License
%along with this program; if not, write to the Free Software
%Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

icatb_defaults;
global DIM_ESTIMATION_OPTS;

if (~exist('files', 'var'))
    files = icatb_selectEntry('typeSelection', 'multiple', 'filter', '*.img', 'title', 'Select functional data files');
end

drawnow;

if (~exist('maskvec', 'var'))
    maskvec = [];
end

if (~exist('precisionType', 'var'))
    precisionType = 'double';
end

dim_est_method = 1;
fwhm = [];

if (~exist('dim_est_opts', 'var') || isempty(dim_est_opts))
    dim_est_opts = DIM_ESTIMATION_OPTS;
end

try
    dim_est_method = dim_est_opts.method;
catch
end


% try
%     dim_est_method = DIM_ESTIMATION_OPTS.method;
% catch
%     if (isfield(DIM_ESTIMATION_OPTS, 'iid_sampling'))
%         if (DIM_ESTIMATION_OPTS.iid_sampling == 1)
%             % iid sampling on
%             dim_est_method = 1;
%         else
%             % iid sampling off
%             dim_est_method = 2;
%         end
%     end
% end

try
    fwhm = dim_est_opts.fwhm;
catch
end

if (dim_est_method == 2)
    if (isempty(fwhm))
        error('Please set fwhm smoothness factor using DIM_ESTIMATION_OPTS.fwhm when DIM_ESTIMATION_OPTS.method = 2');
    end
end

if (ischar(maskvec))
    % If mask is supplied as a file, read data and convert data to
    % logical
    maskvec = icatb_rename_4d_file(maskvec);
    % Load only the first file
    maskvec = icatb_loadData(deblank(maskvec(1, :)));
end

%% Get dimensions of data by loading first time point
if (ischar(files))
    files = icatb_rename_4d_file(files);
    first_file = deblank(files(1, :));
    data = icatb_loadData(first_file);
    xdim = size(data, 1); ydim = size(data, 2); zdim = size(data, 3); tdim = size(files, 1);
    first_file = icatb_parseExtn(first_file);
    if (strcmpi(first_file(end-3:end), '.mat'))
        tdim = size(data, 4);
    end
else
    data = files;
    if (length(size(data)) == 4)
        xdim = size(data, 1); ydim = size(data, 2); zdim = size(data, 3); tdim = size(data, 4);
    else
        if (length(size(maskvec)) == 3)
            xdim = size(maskvec, 1); ydim = size(maskvec, 2); zdim = size(maskvec, 3); tdim = size(data, 2);
            if (size(data, 1) ~= xdim*ydim*zdim)
                error('Mask dimensions (x, y, z) don''t match data dimensions (x, y, z)');
            end
        else
            error('Cannot get the x,y,z dimension information from the mask or the data. Data must be 4D array or mask must be 3D array');
        end
    end
    
end

%% Load mask and check the dimensions of the data
if (isempty(maskvec))
    
    if (length(size(data)) == 4)
        tmp = data(:, :, :, 1);
    elseif ((length(size(data)) == 3) && (numel(data) == xdim*ydim*zdim))
        tmp = data(:);
    elseif ((length(size(data)) == 2) && (numel(data) == xdim*ydim*zdim))
        tmp = data(:);
    else
        tmp = data(:, 1);
    end
    
    
    maskvec = (tmp(:) >= abs(mean(tmp(:))));
    
else
    
    % Get mask vector
    maskvec = (maskvec ~= 0);
    maskvec = maskvec(:);
    
end

if (length(maskvec) ~= xdim*ydim*zdim)
    error('Error:MaskDIM', 'Mask dimensions (%d) doesn''t match the data (%d)', length(maskvec), xdim*ydim*zdim);
end

%% Load data by reading timepoint at a time
if (ischar(files))
    data = icatb_read_data(files, [], find(maskvec == 1), precisionType);
else
    if (strcmpi(precisionType, 'single'))
        data = single(data);
    end
    data = reshape(data, xdim*ydim*zdim, tdim);
    data = data(maskvec, :);
end

verbose = 1;

mdl = [];
aic = [];
kic = [];

if ((dim_est_method == 3) || (dim_est_method == 4))
    % Use ER-FM or ER-AR
    
    if (verbose)
        disp('Using entropy rate estimator ...');
    end
    
    [ER_FM, ER_AR] = orderEstimation_HdsvsEig_R(data');
    
    if (dim_est_method == 3)
        disp('Order estimated by entropy rate based method assumes signal has finite memory length ...');
        comp_est = ER_FM;
    else
        disp('Order estimated by entropy rate based method assumes AR signal ...');
        comp_est = ER_AR;
    end
    
    return;
    
end


if (dim_est_method == 1)
    
    %% Perform variance normalization
    if verbose
        fprintf('\n Performing variance normalization ...');
    end
    
    for n = 1:size(data, 2)
        data(:, n) = detrend(data(:, n), 0) ./ std(data(:, n));
    end
    
    %% (completed)
    
    fprintf('\n P1:');
    [V, EigenValues] = icatb_svd(data, tdim, 'verbose', 0);
    EigenValues = diag(EigenValues);
    EigenValues = EigenValues(end:-1:1);
    dataN = data*V(:, end:-1:1);
    clear V;
    
    %% Select Gaussian principal components for effectively i.i.d. sample estimation
    if verbose
        fprintf('\n Selecting Gaussian principal components ...');
    end
    
    %% Using 12 gaussian components from middle, top and bottom gaussian components to determine the subsampling depth. Final subsampling depth is determined using median
    kurtv1 = kurtn(dataN);
    kurtv1(EigenValues > mean(EigenValues)) = 1000;
    idx_gauss = find(((kurtv1(:) < 0.3) .* (EigenValues > eps)) == 1);
    idx = idx_gauss(:)';
    dfs = length(find(EigenValues > eps)); % degrees of freedom
    minTp = 12;
    
    if (length(idx) >= minTp)
        middle = round(length(idx)/2) + 1;
        idx = [idx(1:4), idx(middle - 1:middle + 2), idx(end-3:end)];
    elseif isempty(idx)
        minTp = min([minTp, dfs]);
        idx = (dfs - minTp + 1:dfs);
    end
    
    idx = unique(idx);
    
    %% Estimate the subsampling depth for effectively i.i.d. samples
    if verbose
        fprintf('\n\n Estimate effectively i.i.d. samples: ');
    end
    
    mask_ND = reshape(maskvec, xdim, ydim, zdim);
    ms = length(idx);
    s = zeros(ms,1);
    for i = 1:ms
        x_single = zeros(xdim*ydim*zdim, 1);
        x_single(maskvec) = dataN(:, idx(i));
        x_single = reshape(x_single, xdim, ydim, zdim);
        s(i) = est_indp_sp(x_single);
        if (i > 6)
            tmpS = s(1:i);
            tmpSMedian = round(median(tmpS));
            if (length(find(tmpS == tmpSMedian)) > 6)
                s = tmpS;
                break;
            end
        end
        dim_n = prod(size(size(x_single)));
        clear x_single;
    end
    clear dataN;
    fprintf('\n');
    s1 = round(median(s));
    if floor((sum(maskvec)/(1*tdim))^(1/dim_n)) < s1
        s1 = floor((sum(maskvec)/(1*tdim))^(1/dim_n));
    end
    N = round(sum(maskvec)/s1^dim_n);
    %% (completed)
    
    %% Use the subsampled dataset to calculate eigen values
    if verbose
        fprintf('\n Perform EVD on the effectively i.i.d. samples ...');
    end
    if s1 ~= 1
        mask_s = subsampling(mask_ND, s1, [1,1,1]);
        mask_s_1d = reshape(mask_s, 1, prod(size(mask_s)));
        dat = zeros(length(find(mask_s_1d == 1)), tdim);
        for i = 1:tdim
            x_single = zeros(xdim*ydim*zdim, 1);
            x_single(maskvec) = data(:, i);
            x_single = reshape(x_single, xdim, ydim, zdim);
            dat0 = subsampling(x_single, s1, [1, 1, 1]);
            dat0 = reshape(dat0, prod(size(dat0)), 1);
            dat(:, i) = dat0(mask_s_1d);
            clear x_single;
        end
        clear data;
        %% Perform variance normalization
        for n = 1:size(dat, 2)
            dat(:, n) = detrend(dat(:, n), 0) ./ std(dat(:, n));
        end
        %% (completed)
        fprintf('\n P2:');
        [V, EigenValues] = icatb_svd(dat, tdim, 'verbose', 0);
        EigenValues = diag(EigenValues);
        EigenValues = EigenValues(end:-1:1);
        %lam = EigenValues;
    end
    
    lam = EigenValues;
    clear dat;
    
    if verbose
        fprintf('\n Effective number of i.i.d. samples: %d \n',N);
    end
    
    %% Make eigen spectrum adjustment
    if verbose
        fprintf('\n Perform eigen spectrum adjustment ...');
    end
    lam = eigensp_adj(lam, N, length(lam));
    %% (completed)
    
    if sum(imag(lam))
        error('Invalid eigen value found for the subsampled data.');
    end
    
    
else
    %% No iid sampling
    if verbose
        fprintf('\n computing eigen values ...\n');
    end
    
    [V, EigenValues] = icatb_svd(data, tdim, 'verbose', 0);
    EigenValues = diag(EigenValues);
    EigenValues = EigenValues(end:-1:1);
    lam = EigenValues;
    N = ceil(size(data, 1) / prod(fwhm));
    
end

%% Correction on the ill-conditioned results (when tdim is large, some least significant eigenvalues become small negative numbers)
lam(real(lam) <= eps) = min(lam(real(lam)>= eps));

if verbose
    fprintf('\n Estimating the dimension ...');
end
p = tdim;
aic = zeros(1, p - 1);
kic = zeros(1, p - 1);
mdl = zeros(1, p - 1);
for k = 1:p-1
    LH = log(prod( lam(k+1:end).^(1/(p-k)) )/mean(lam(k+1:end)));
    mlh = 0.5*N*(p-k)*LH;
    df = 1 + 0.5*k*(2*p-k+1);
    aic(k) =  -2*mlh + 2*df;
    kic(k) =  -2*mlh + 3*df;
    mdl(k) =  -mlh + 0.5*df*log(N);
end

% Find the first local minimum of each ITC
itc = zeros(3, length(mdl));
itc(1,:) = aic;
itc(2,:) = kic;
itc(3,:) = mdl;

%% Estimated components using MDL
dlap = squeeze(itc(3, 2:end) - itc(3, 1:end-1));
a = find(dlap > 0);
if isempty(a)
    comp_est = length(squeeze(itc(3, :)));
else
    comp_est = a(1);
end

fprintf('\n');
disp(['Estimated components is found out to be ', num2str(comp_est)]);
fprintf('\n');

function out = subsampling(x,s,x0)
% Subsampling the data evenly with space 's'

n = size(x);

if max(size(n)) == 2 & min(n) == 1  % 1D
    out = x([x0(1):s:max(n)]);
    
elseif max(size(n)) == 2 & min(n) ~= 1  % 2D
    out = x([x0(1):s:n(1)],[x0(2):s:n(2)]);
    
elseif max(size(n)) == 3 & min(n) ~= 1  % 3D
    out = x([x0(1):s:n(1)],[x0(2):s:n(2)],[x0(3):s:n(3)]);
    
else
    error('Unrecognized matrix dimension!(subsampling)');
    
end


function [s, entrate_m] = est_indp_sp(x)
% estimate the effective number of independent samples based on the maximum
% entropy rate principle of stationary random process


dimv = size(x);
s0 = 0;
fprintf('\n');
for j = 1:min(dimv)-1
    x_sb = subsampling(x,j,[1,1,1]);
    if j == 1
        fprintf('\n Estimating the entropy rate of the Gaussian component with subsampling depth %d,',j);
    else
        fprintf(' %d,',j);
    end
    entrate_m = entrate_sp(x_sb,1);
    
    ent_ref = 1.41;
    if entrate_m > ent_ref
        s0 = j; break;
    end
end
fprintf(' Done;');
if s0 == 0
    error('Ill conditioned data, can not estimate independent samples.(est_indp_sp)');
else
    s = s0;
end


function out = entrate_sp(x, sm_window)
% Calculate the entropy rate of a stationary Gaussian random process using
% spectrum estimation with smoothing window

n = size(x);

% Normalize x_sb to be unit variance
x_std = std(reshape(x,prod(n),1));
if x_std < 1e-10; x_std = 1e-10; end;
x = x/x_std;

if( sm_window == 1)
    
    M = ceil(n/10);
    
    if(max(size(n) >= 3))
        parzen_w_3 = zeros(2*n(3)-1,1);
        parzen_w_3(n(3)-M(3):n(3)+M(3)) = parzen_win(2*M(3)+1);
    end
    
    if(max(size(n) >= 2))
        parzen_w_2 = zeros(2*n(2)-1,1);
        parzen_w_2(n(2)-M(2):n(2)+M(2)) = parzen_win(2*M(2)+1);
    end
    
    if(max(size(n) >= 1))
        parzen_w_1 = zeros(2*n(1)-1,1);
        parzen_w_1(n(1)-M(1):n(1)+M(1)) = parzen_win(2*M(1)+1);
    end
    
end

if max(size(n)) == 2 & min(n) == 1  % 1D
    xc = xcorr(x,'unbiased');
    xc = xc.*parzen_w;
    xf = fftshift(fft(xc));
    
elseif max(size(n)) == 2 & min(n) ~= 1  % 2D
    xc = xcorr2(x); % default option: computes raw correlations with NO
    % normalization -- Matlab help on xcorr
    
    % Bias correction
    v1 = [1:n(1),n(1)-1:-1:1];
    v2 = [1:n(2),n(2)-1:-1:1];
    
    vd = v1'*v2;
    xc = xc./vd;
    parzen_window_2D = parzen_w_1*parzen_w_2';
    xc = xc.*parzen_window_2D;
    xf = fftshift(fft2(xc));
    
elseif max(size(n)) == 3 & min(n) ~= 1  % 3D
    xc = zeros(2*n-1);
    for m3 = 0:n(3)-1
        temp = zeros(2*n(1:2)-1);
        for k = 1:n(3)-m3
            temp = temp + xcorr2(x(:,:,k+m3),x(:,:,k)); % default option:
            % computes raw correlations with NO normalization
            % -- Matlab help on xcorr
        end
        xc(:,:,n(3)-m3) = temp;
        xc(:,:,n(3)+m3) = temp;
    end
    
    % Bias correction
    v1 = [1:n(1),n(1)-1:-1:1];
    v2 = [1:n(2),n(2)-1:-1:1];
    v3 = [n(3):-1:1];
    
    vd = v1'*v2;
    vcu = zeros(2*n-1);
    for m3 = 0:n(3)-1
        vcu(:,:,n(3)-m3) = vd*v3(m3+1);
        vcu(:,:,n(3)+m3) = vd*v3(m3+1);
    end
    
    xc = xc./vcu;
    
    parzen_window_2D = parzen_w_1*parzen_w_2';
    for m3 = 0:n(3)-1
        parzen_window_3D(:,:,n(3)-m3) = parzen_window_2D*parzen_w_3(n(3)-m3);
        parzen_window_3D(:,:,n(3)+m3) = parzen_window_2D*parzen_w_3(n(3)+m3);
    end
    xc = xc.*parzen_window_3D;
    
    xf = fftshift(fftn(xc));
    
else
    error('Unrecognized matrix dimension.');
    
end

xf = abs(xf);
xf(xf<1e-4) = 1e-4;
out = 0.5*log(2*pi*exp(1)) + sumN(log(abs((xf))))/2/sumN(abs(xf));


function w = parzen_win(n)
% PARZENWIN Parzen window.
%   PARZENWIN(N) returns the N-point Parzen (de la Valle-Poussin) window in
%   a column vector.

% Check for valid window length (i.e., n < 0)
[n,w,trivialwin] = checkOrder(n);
if trivialwin, return, end;

% Index vectors
k = -(n-1)/2:(n-1)/2;
k1 = k(k<-(n-1)/4);
k2 = k(abs(k)<=(n-1)/4);

% Equation 37 of [1]: window defined in three sections
w1 = 2 * (1-abs(k1)/(n/2)).^3;
w2 = 1 - 6*(abs(k2)/(n/2)).^2 + 6*(abs(k2)/(n/2)).^3;
w = [w1 w2 w1(end:-1:1)]';


function [n_out, w, trivalwin] = checkOrder(n_in)
% CHECKORDER Checks the order passed to the window functions.

w = [];
trivalwin = 0;

% Special case of negative orders:
if n_in < 0,
    error('Order cannot be less than zero.');
end

% Check if order is already an integer or empty
% If not, round to nearest integer.
if isempty(n_in) | n_in == floor(n_in),
    n_out = n_in;
else
    n_out = round(n_in);
    warning('Rounding order to nearest integer.');
end

% Special cases:
if isempty(n_out) | n_out == 0,
    w = zeros(0,1);               % Empty matrix: 0-by-1
    trivalwin = 1;
elseif n_out == 1,
    w = 1;
    trivalwin = 1;
end


function lam_adj = eigensp_adj(lam,n,p)
% Eigen spectrum adjustment for EVD on finite samples

r = p/n;

bp = (1+sqrt(r))^2;
bm = (1-sqrt(r))^2;

vv = [bm:(bp-bm)/(5*p-1):bp];

gv = 1./(2*pi*r*vv).*sqrt((vv-bm).*(bp-vv));

gvd = zeros(size(gv));
for i = 1:length(gv)
    gvd(i) = sum(gv(1:i));
end
gvd = gvd/max(gvd);

lam_emp = zeros(size(lam));
for i = 1:p
    i_norm = i/p;
    [minv,minx] = min(abs(i_norm-gvd));
    lam_emp(i) = vv(minx);
end

lam_emp = rot90(lam_emp, 2);

lam_adj = lam./lam_emp;


function [sum_dat] = sumN(dat)
% sum of the all elements of the dat matrix

sum_dat = sum(dat(:));


function c = xcorr2(a,b)
% two dimensional cross correlation

if nargin == 1
    b = a;
end

c = conv2(a, rot90(conj(b),2));

function kurt = kurtn(x)
% Normalized kurtosis funtion so that for a Gaussian r.v. the kurtn(g) = 0
% Input: x 1:N vector

kurt = zeros(1, size(x, 2));

for i = 1:size(x, 2)
    a = detrend(x(:, i), 0);
    a = a/std(a);
    kurt(i) = mean(a.^4) - 3;
end



function [ER_FM, ER_AR] = orderEstimation_HdsvsEig_R( data1 )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ENTROPY-RATE BASED ORDER SELECTION (ER_FM & ER_AR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input:
% data1:       mixtures X
% Outputs:
% ER_FM:       Order estimated by entropy rate based method that assumes signal has finite memory length
% ER_AR:       Order estimated by entropy rate based method that assumes AR signal

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% References:
%
% ER_FM & ER_AR
% G.-S. Fu, M. Anderson, and T. Adali, Likelihood estimators for dependent samples
% and their application to order detection, in Signal Process, IEEE Trans.
% on, vol. 62, no. 16, pp. 4237-4244, Aug. 2014.
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program by Geng-Shen Fu
% Please contact me at fugengs1@umbc.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data1 = data1 - mean(data1,2)*ones(1,size(data1,2));    % remove the mean
[N T]=size(data1);

% estimate the downsampling depth and order detection
[U1,S1,V1] = svd(data1,'econ');
clear V1;
data1 = U1'*data1;

[ACS DSD dsdAll] = Entropy_rate_estimator(data1);

ER_FM_DS=zeros(1,max(DSD,100));
for i=1:DSD
    for n=1:N
        H(n)=Entropy_rate_single_estimator(ACS(n,1:i:DSD));
    end;
    if sum(H==Inf)>0
        continue;
    end
    [H index_H] = sort(H, 'descend');
    ER_FM_DS(i)=order_estimate_entropy(N, T/i, H);
end
ER_FM_DS1=ER_FM_DS(ER_FM_DS~=0);
ER_FM=round(median(ER_FM_DS1));

ER_AR_DS=zeros(1,max(DSD,100));
for i=1:DSD
    for n=1:N
        H(n)=Entropy_rate_single_estimator_AR(data1(n,1:i:end));
    end;
    if sum(H==Inf)>0
        continue;
    end
    [H index_H] = sort(H, 'descend');
    ER_AR_DS(i)=order_estimate_entropy(N, T/i, H);
end
ER_AR_DS1=ER_AR_DS(ER_AR_DS~=0);
ER_AR=round(median(ER_AR_DS1));

maxp0 = floor(2*size(data1,2)/size(data1,1))/2;
dsd = size(data1,1)/2;    % the initial downsampling depth is half the number of sensors
existing_orders = [];   % save all estimated orders
while 1
    % estimate orders
    E_JDS = order_estimate_complex_c(size(data1,1), size(data1,2), dsd, diag(S1));
    existing_orders = [existing_orders; E_JDS];
    if sum(existing_orders==E_JDS)>=5  % return when the order detection converges
        dsd = ceil(dsd);
        break;
    end
    
    % re-estimate downsampling depth
    dsd = 0;
    cnt1 = 0;
    cnt2 = 0;   % counting the MA orders reaching the upper bound
    ma_order_list=[];
    for k=1:E_JDS
        ma_order1=DownSampOrderEstimationR(data1(k,:));
        ma_order_list=[ma_order_list;ma_order1];
    end
    dsd=median(ma_order_list);
end

function p0 = order_estimate_complex_c( N, T0, dsd, S1 )
%  order estimation of real-valued multivariate process
%  N: number of time points
%  T0: original sample size
%  dsd: downsampling depth
%  S1: singular values
%  p0: estimated order
% Program by Geng-Shen Fu: fugengs1@umbc.edu
T = T0/dsd;     % effective sample size
lambda = S1.^2/T0;   % eigenvalues

DL = zeros(N,1);
for p = 0 : N-1
    
    nf = 1+p*N-0.5*p*(p-1); % number of degrees of freedom
    DL(p+1) = DL(p+1) + 0.5*nf*log(T)/T;
    
    for r = 1 : p
        DL(p+1) = DL(p+1) + 0.5*log(lambda(r));
    end
    
    sigma2v = sum(lambda(p+1:N))/(N-p);
    DL(p+1) = DL(p+1) + 0.5*(N-p)*log(sigma2v);
    
end
[dl0, p0] = min(real(DL));
p0 = p0-1;
%figure; plot([0:N-1], real(DL),'--o')



function [ACS dsd ma_order_list] = Entropy_rate_estimator( data )
N = size(data,1);
T = size(data, 2);
ma_order_list=[];

for i = 1:N
    dsd=DownSampOrderEstimationR(data(i,:));
    ma_order_list = [ma_order_list;dsd];
end;
% ma_order_list=[40 35 30 25 10 10 10 10 10 10]';

dsd=round(mean(ma_order_list));
dsd_max=max(ma_order_list);

% ma_order_list(ma_order_list<dsd)=dsd;
% dsd=105;
for i=1:N
    data(i,:)=data(i,:)-mean(data(i,:));
    acsTmp=xcorr(data(i,:),dsd_max-1,'unbiased');
    ACS(i,:)=acsTmp(dsd_max:end);
end;


function [H1 E1 C1 C2 ma_order_list] = Entropy_estimate_complex( data )
N = size(data,1);
T = size(data, 2);
ma_order_list=[];

for i = 1:N
    dsd=DownSampOrderEstimationR(data(i,:));
    ma_order_list = [ma_order_list;dsd];
end;
% ma_order_list=[40 35 30 25 10 10 10 10 10 10]';

dsd_max = max(ma_order_list);
for i = 1:N
    [H1(i) C1(:,:,i) C2(:,:,i)] = Entropy_single_estimate_complex(data(i,:), ma_order_list(i), dsd_max);
    E1(i) = var(data(i,1:dsd:T), 1);
end;

function [H] = Entropy_rate_single_estimator( acs )
% Program by Geng-Shen Fu: fugengs1@umbc.edu
dsd=size(acs,2);
C1=toeplitz(acs);
if dsd~=1
    C2=toeplitz(acs(1:end-1));
else
    C2=1;
end;
acs1=acs(end:-1:2);
if det(C2)<1e-18
    H=Inf;
    return;
end
H= acs(1)-acs1*inv(C2)*acs1';
% H= (abs(det(C1))/abs(det(C2)));


function [H] = Entropy_rate_single_estimator_AR( data )
% Program by Geng-Shen Fu: fugengs1@umbc.edu
%AR Order Selection with Partial Autocorrelation Sequence
T=size(data,2);
[arcoefs,E,K] = aryule(data,15);

AR_Order=find(abs(K)<1.96/sqrt(T),1)-1;
if prod(size(AR_Order))==0
    AR_Order=15;
end
[arcoefs] = aryule(data,AR_Order);
Data_White = filter(arcoefs, 1, data);
H=var(Data_White);

% figure;
% pacf = -K;
% lag = 1:15;
% stem(lag,pacf,'markerfacecolor',[0 0 1]);
% xlabel('Lag'); ylabel('Partial Autocorrelation');
% set(gca,'xtick',1:1:15)
% lconf = -1.96/sqrt(T)*ones(length(lag),1);
% uconf = 1.96/sqrt(T)*ones(length(lag),1);
% hold on;
% line(lag,lconf,'color',[1 0 0]);
% line(lag,uconf,'color',[1 0 0]);
% i=1;



function [H1 C1 C2] = Entropy_single_estimate_complex( data, dsd, dsd_max )
% Program by Geng-Shen Fu: fugengs1@umbc.edu
data=data-mean(data);
acs=xcorr(data,dsd-1,'unbiased');
C1=toeplitz(acs(dsd:end));
if dsd~=1
    C2=toeplitz(acs(dsd:end-1));
else
    C2=1;
end;
H1= (abs(det(C1))/abs(det(C2)));

acs=xcorr(data,dsd_max-1,'unbiased');
C1=toeplitz(acs(dsd_max:end));
if dsd_max~=1
    C2=toeplitz(acs(dsd_max:end-1));
else
    C2=1;
end;

function [p0] = order_estimate_entropy( N, T0, H1 )
%  order estimation of real-valued multivariate process
%  N: number of time points
%  T0: original sample size
%  dsd: downsampling depth
%  S1: singular values
%  p0: estimated order
% Program by Geng-Shen Fu: fugengs1@umbc.edu
% T = T0/dsd;     % effective sample size
T = T0;     % effective sample size
lambda = H1;   % eigenvalues

DL = zeros(N,1);
for p = 0 : N-1
    nf = 1+p*N-0.5*p*(p-1); % number of degrees of freedom
    
    for r = 1 : p
        DL(p+1) = DL(p+1) + 0.5*log(lambda(r));
    end
    
    sigma2v = sum(lambda(p+1:N))/(N-p);
    DL(p+1) = DL(p+1) + 0.5*(N-p)*log(sigma2v);
    
    DL(p+1) = DL(p+1) + 0.5*nf*log(T)/T;
end
[dl0, p0] = min(real(DL));
p0 = p0-1;
% DL'
%figure; plot([0:N-1], real(DL),'--o')

function [p0] = order_estimate_covariance( N, T, C1, C2, dsd )
%  order estimation of real-valued multivariate process
%  N: number of time points
%  T0: original sample size
%  dsd: downsampling depth
%  S1: singular values
%  p0: estimated order
% Program by Geng-Shen Fu: fugengs1@umbc.edu

DL = zeros(N,1);
for p = 0 : N-1
    nf = p*(2*N-p-1)/2+(p+1)*dsd; % number of degrees of freedom
    
    for r = 1 : p
        DL(p+1) = DL(p+1) + 0.5*log((abs(det(C1(:,:,r)))/abs(det(C2(:,:,r)))));
    end
    
    C1Bar = sum(C1(:,:,p+1:N),3)/(N-p);
    C2Bar = sum(C2(:,:,p+1:N),3)/(N-p);
    DL(p+1) = DL(p+1) + 0.5*(N-p)*log((abs(det(C1Bar))/abs(det(C2Bar))));
    
    DL(p+1) = DL(p+1) + 0.5*nf*log(T)/T;
end
[dl0, p0] = min(real(DL));
p0 = p0-1;
% DL'
%figure; plot([0:N-1], real(DL),'--o')


%For debug
function DisplayACS(data)
lag=round(prod(size(data))/20);
acs=abs(xcorr(data, lag, 'unbiased'));
figure;
plot(1:lag+1,acs(lag+1:end));




function p0 = DownSampOrderEstimationR(data)
% return the estimated order of moving average process x
x = data - mean(data);
x = x/std(x);
p0 = 1;

while 1
    if isWhiteR( x(ceil(rand*p0):p0:end) ) % ceil(rand*p0) is a random offset
        return;
    else
        p0=p0+1;
        if p0>inf
            p0=inf;
            return;
        end
    end
end


function s = isWhiteR(x)
% hypothesis: white process vs colored process
% if accept white, return 1; otherwise, return 0
% Program by Xi-Lin Li: lixilin@umbc.edu
T=length(x);

% DL of white case
p=0;
dl0=log(mean(abs(x).^2));

% DL of colored case with AR order 1
p=1;
x1 = [0,x(1:end-1)];
a = x*x1'/(x1*x1');
n = x-a*x1;
dl1 = log(mean(abs(n).^2)) + p*log(T)/T/2;

s=(dl0<dl1);
