function [comp_est, mdl, aic, kic] = icatb_estimate_dimension(files, maskvec, precisionType)
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
%
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

iid_sampling = 1;
fwhm = [];
try
    iid_sampling = DIM_ESTIMATION_OPTS.iid_sampling;
catch
end

try
    fwhm = DIM_ESTIMATION_OPTS.fwhm;
catch
end

if (~iid_sampling)
    if (isempty(fwhm))
        error('Please set fwhm smoothness factor when i.i.d sampling is turned off in DIM_ESTIMATION_OPTS');
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


if (iid_sampling)
    
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