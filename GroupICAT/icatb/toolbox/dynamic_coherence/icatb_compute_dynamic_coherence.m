function [Corrected_Out_Distribution, Null_Distribution, mask_v_inds, coin] = icatb_compute_dynamic_coherence(TC, settings_file)
%% Dynamic coherence
%
% Inputs:
% 1. TC - Timecourses matrix
% 2. settings_file - Settings MAT file
%
% Outputs:
% 1. Corrected_Out_Distribution -
% 2. mask_v_inds -
% 3. coin -
%
%


useWavToolbox = 1;
try
    wmaxlev(2^10, 'db1');
catch
    useWavToolbox = 0;
end

load(settings_file);


TC_ICN = TC;

for n = 1:size(TC_ICN, 1)
    tmp = squeeze(TC_ICN(n, :, :));
    tmp = icatb_zscore(tmp);
%     tmp = zscores(tmp);
%     tmp = tmp';
    TC_ICN(n, :, :) = tmp;
end

L = size(TC_ICN, 2);

Real_Time = 0:L-1;
Real_Time = dt.*Real_Time;

Psi = computeWavelet(squeeze(TC_ICN(1,:,1)), 0.5/(Real_Time(2)-Real_Time(1)), fscales);


TC_ICN2 = zeros(size(TC_ICN, 1), size(TC_ICN, 3), length(fscales), size(TC_ICN, 2));
for n = 1:size(TC_ICN2, 1)
    for nComp = 1:size(TC_ICN2, 2)
        TC_ICN2(n, nComp, :, :) = my_CWT(squeeze(TC_ICN(n, :, nComp)), fscales, 0.5/(Real_Time(2)-Real_Time(1)), Psi);
    end
end


coin = my_COI(L,fscales,dt,useWavToolbox);
coin = 1 - coin;

[Null_Distribution, dd] = Null_Wavelet_SlidingWindow_TimeIndependent(TC_ICN, Real_Time, fscales, my_ntr,my_ntr_2, [0, 1], [], coin, Psi, TC_ICN2);

[Out_Distribution, mask_v_inds] = Wavelet_SlidingWindow(TC_ICN, Real_Time, fscales, my_ntr,my_ntr_2, [0, 1], [], coin, Psi);

Null_Distribution = repmat(Null_Distribution, [1, 1, 1, size(coin,2) ]);


Null_Distribution = reshape(Null_Distribution, size(Null_Distribution,1), size(Null_Distribution,2), []);
Mean_Null_Dist = Null_Distribution(1, :, :);
Std_Null_Dist = Null_Distribution(2, :, :);


Corrected_Out_Distribution = Out_Distribution - repmat(Mean_Null_Dist, [size(Out_Distribution,1), 1, 1]);

Std_Null_Dist = repmat(Std_Null_Dist, [size(Out_Distribution,1), 1, 1]);
Corrected_Out_Distribution = real(Corrected_Out_Distribution) ./ real(Std_Null_Dist) + ...
    1i*imag(Corrected_Out_Distribution) ./ imag(Std_Null_Dist);



function [Out_Distribution mask_v_inds incoi_mask] = Wavelet_SlidingWindow(TC_ICN, Real_Time, scales, my_ntr,my_ntr_2, freq_smoothing, Null_Distribution, incoi_mask, Psi)
%% Sliding Windows: win_sz
%     nC = Num_of_Clusters;
% CovarianceMatrices = atanh(0.999.*CovarianceMatrices);
%

% Surrogate TC_ICN
% TC_ICN_v = permute(TC_ICN, [1 3 2]);
% TC_ICN_v = reshape(TC_ICN_v, size(TC_ICN_v,1)*size(TC_ICN_v,2), []);
% % for i=1:size(TC_ICN_v,1)
% %     TC_ICN_v(i,:) = Surrogate_FourierTransform(TC_ICN_v(i,:));
% % end
% TC_ICN_v = Surrogate_FourierTransform(TC_ICN_v')';
%
% TC_ICN_v = reshape(TC_ICN_v, size(TC_ICN,1), size(TC_ICN,3), []);
% TC_ICN = permute(TC_ICN_v, [1 3 2]);

alaki_mask = ones(size(TC_ICN,3));
mask_inds = alaki_mask - triu(alaki_mask);
mask_v_inds = find(mask_inds > 0);

L = length(Real_Time);


incoi_mask_v = find(incoi_mask > 0);

dt = Real_Time(2)-Real_Time(1);


%Out_Distribution = zeros(size(TC_ICN,1), length(mask_v_inds), length(incoi_mask(:))) + 1i*zeros(size(TC_ICN,1), length(mask_v_inds), length(incoi_mask(:)));
%Out_Distribution = rand(size(TC_ICN,1), length(mask_v_inds), length(incoi_mask_v)) + 1i*rand(size(TC_ICN,1), length(mask_v_inds), length(incoi_mask_v));
%Out_Distribution = rand(nchoosek(size(TC_ICN,1),2), length(mask_v_inds), length(incoi_mask_v));

Out_Distribution = zeros(nchoosek(size(TC_ICN,1),2), length(mask_v_inds), length(incoi_mask(:)))  + ...
    1i*zeros(nchoosek(size(TC_ICN,1),2), length(mask_v_inds), length(incoi_mask(:)));

current_subject_step = 1;
for subject_num_one = 1:size(TC_ICN,1)%rand_inds
    for subject_num_two = subject_num_one+1:size(TC_ICN,1)%rand_inds
        fprintf('Processing Subject # %d out of %d\n',current_subject_step, size(Out_Distribution,1));
        
        current_comp_step = 1;
        for comp_one = 1:size(TC_ICN,3)
            %fprintf('Processing Subject # %d Comp %d\n',subject_num, comp_one);
            for comp_two = comp_one+1:size(TC_ICN,3)
                
                %                 first_timeseries = [Real_Time(:), zscore(squeeze(TC_ICN(subject_num_one,:,comp_one)))'];
                %                 second_timeseries = [Real_Time(:), zscore(squeeze(TC_ICN(subject_num_two,:,comp_two)))'];
                True_Signal(1,:) = squeeze(TC_ICN(subject_num_one,:,comp_one));
                True_Signal(2,:) = squeeze(TC_ICN(subject_num_two,:,comp_two));
                
                %True_Signal = zscores(True_Signal);
                
                
                tmp_coeffs = tfcohnewwavelet3(True_Signal(1,:), True_Signal(2,:),scales,my_ntr,my_ntr_2,freq_smoothing(1),freq_smoothing(2),dt,incoi_mask, Null_Distribution, Psi);
                
                
                %
                Out_Distribution (current_subject_step,current_comp_step,incoi_mask_v) = tmp_coeffs(incoi_mask_v);%tmp_coeffs(incoi_mask_v);%sum(abs(tmp_masked),2)./sum(incoi_mask,2);
                current_comp_step = current_comp_step + 1;
                
            end
        end
        current_subject_step = current_subject_step + 1;
    end
end



function [Out_Distribution mask_v_inds] = Null_Wavelet_SlidingWindow_TimeIndependent(TC_ICN, Real_Time, scales, my_ntr,my_ntr_2, freq_smoothing, Null_Distribution, incoi_mask, Psi, TC_ICN2)
%% Sliding Windows: win_sz
%     nC = Num_of_Clusters;
% CovarianceMatrices = atanh(0.999.*CovarianceMatrices);
%

% Surrogate TC_ICN
% TC_ICN_v = permute(TC_ICN, [1 3 2]);
% TC_ICN_v = reshape(TC_ICN_v, size(TC_ICN_v,1)*size(TC_ICN_v,2), []);
% % for i=1:size(TC_ICN_v,1)
% %     TC_ICN_v(i,:) = Surrogate_FourierTransform(TC_ICN_v(i,:));
% % end
% TC_ICN_v = Surrogate_FourierTransform(TC_ICN_v')';
%
% TC_ICN_v = reshape(TC_ICN_v, size(TC_ICN,1), size(TC_ICN,3), []);
% TC_ICN = permute(TC_ICN_v, [1 3 2]);

alaki_mask = ones(size(TC_ICN,3));
mask_inds = alaki_mask - triu(alaki_mask);
mask_v_inds = find(mask_inds > 0);

L = length(Real_Time);


incoi_mask_v = 1:size(incoi_mask,1);%find(incoi_mask > 0);


time_samples_sz = min(sum(incoi_mask,2), round(40*size(incoi_mask,2)/100));

dt = Real_Time(2)-Real_Time(1);


%Out_Distribution = zeros(size(TC_ICN,1), length(mask_v_inds), length(incoi_mask(:))) + 1i*zeros(size(TC_ICN,1), length(mask_v_inds), length(incoi_mask(:)));
%Out_Distribution = rand(size(TC_ICN,1), length(mask_v_inds), length(incoi_mask_v)) + 1i*rand(size(TC_ICN,1), length(mask_v_inds), length(incoi_mask_v));
Out_Distribution = zeros(2, length(mask_v_inds), length(incoi_mask_v))  + ...
    1i*zeros(2, length(mask_v_inds), length(incoi_mask_v));
Tmp_Out_Distribution = zeros(length(incoi_mask_v),nchoosek(size(TC_ICN,1),2),max(time_samples_sz))  + ...
    1i*zeros(length(incoi_mask_v),nchoosek(size(TC_ICN,1),2),max(time_samples_sz));
current_cmppair_step = 1;
for comp_one = 1:size(TC_ICN,3)
    %fprintf('Processing Subject # %d Comp %d\n',subject_num, comp_one);
    for comp_two = comp_one+1:size(TC_ICN,3)
        
        
        fprintf('Processing Component pair # %d out of %d\n',current_cmppair_step, size(Out_Distribution,2));
        current_sbjpair_step = 1;
        for subject_num_one = 1:size(TC_ICN,1)%rand_inds
            for subject_num_two = subject_num_one+1:size(TC_ICN,1)%rand_inds
                
                
                %                 first_timeseries = [Real_Time(:), zscore(squeeze(TC_ICN(subject_num_one,:,comp_one)))'];
                %                 second_timeseries = [Real_Time(:), zscore(squeeze(TC_ICN(subject_num_two,:,comp_two)))'];
                %True_Signal(1,:) = squeeze(TC_ICN(subject_num_one,:,comp_one));
                %True_Signal(2,:) = squeeze(TC_ICN(subject_num_two,:,comp_two));
                
                %True_Signal = zscores(True_Signal);
                
                
                tmp_coeffs = tfcohnewwavelet3(squeeze(TC_ICN2(subject_num_one, comp_one, :, :)), squeeze(TC_ICN2(subject_num_two, comp_two, :, :)), scales,my_ntr,my_ntr_2,freq_smoothing(1),freq_smoothing(2),dt,incoi_mask, Null_Distribution, Psi);
                
                for freq_ind =1:size(tmp_coeffs,1)
                    %
                    tmp_validinds = find(incoi_mask(freq_ind,:));
                    rnd_timepoints = randi([1 length(tmp_validinds)], [time_samples_sz(freq_ind), 1]);
                    Tmp_Out_Distribution (freq_ind,current_sbjpair_step,1:length(rnd_timepoints)) = (tmp_coeffs(freq_ind,tmp_validinds(rnd_timepoints)));%tmp_coeffs(incoi_mask_v);%sum(abs(tmp_masked),2)./sum(incoi_mask,2);
                end
                
                current_sbjpair_step = current_sbjpair_step + 1;
            end
        end
        
        for freq_ind = 1:size(Tmp_Out_Distribution,1)
            current_samples = squeeze(Tmp_Out_Distribution(freq_ind,:,1:time_samples_sz(freq_ind)));
            tmp_mean = mean(current_samples(:));
            tmp_std = std(real(current_samples(:))) + ...
                1i*std(imag(current_samples(:)));
            Out_Distribution(1,current_cmppair_step,freq_ind) = tmp_mean;
            Out_Distribution(2,current_cmppair_step,freq_ind) = tmp_std;
        end
        
        current_cmppair_step = current_cmppair_step + 1;
    end
end


function [coi] = my_COI(L,center_freqs, dT, useWavToolbox)
Fs = 0.5/dT;
L = L+mod(L,2);


coi_radius = round(1./(dT*center_freqs));
coi_radius = coi_radius - mod(coi_radius,2);
coi_radius = coi_radius/2;
coi = ones(length(center_freqs), L);




for curr_step=1:length(center_freqs)
    wname = sprintf('cmor%d-%d',round(L/10),center_freqs(curr_step));
    %[Psi,time_points] = cmorwavf(-dT*L/2,dT*L/2,L,round(L/4),center_freqs(curr_step));
    
    if (useWavToolbox)
        [~,l_inds] = find(conofinf(wname,1,round(dT*L),0));
    else
        [~,l_inds] = find(getcoi(1,round(dT*L),0));
    end
    
    if (useWavToolbox)
        [~,h_inds] = find(conofinf(wname,1,round(dT*L),round(dT*L)+1));
    else
        [~,h_inds] = find(getcoi(1,round(dT*L),round(dT*L)+1));
    end
    
    coi(curr_step,max(round(l_inds/dT))+coi_radius(curr_step)+1:min(round(h_inds/dT))-1-coi_radius(curr_step)) = 0;
    
    
end


function varargout = tfcohnewwavelet3(x_new,y_new,fscales,time_winrad1, time_winrad2, scale_winrad1, scale_winrad2, dt, coin, Null_Distribution, Psi)
% TFCOHF3 Time-frequency coherency
% Estimates the complex coherency coefficients using Fourier decomposition
% of vector X and vector Y. The cross and auto spectra are smoothed with
% non-identical smoothing windows (sm_win1 and sm_win2).
%
% Time-frequency coherency is computed by smoothing the cross- and autospectra
% using a smoothing kernel specficied by SM_WIN. The cross- and autospectra
% are estimated using Welch's periodogram method. Signals are dived into overlapping
% sections, each of which is detrended and windowed by the SPEC_WIN parameter,
% then zero padded to length NFFT. TSTEP defines the number of samples the
% window is slided forward. The spectra X, Y, and XY are estimated for each
% segment.Spectral coefficients are then smoothed for the estimation of
% time-frequency coherency using identical smoothing windows.
%
% ARGUMENTS:
%           x           --  signal 1 (vector)
%           y           --  signal 2 (vector)

%                           than n points

%           sm_win1     --  length of window used for smoothing the cross
%                           spectrum (Gauss window)
%           sm_win2     --  length of window used for smoothing the auto
%                           spectra (Gauss window)

%           dt          --  sample rate
%
% If length(sm_win)==1 it specifies the length of the smoothing window in
% seconds. If length(sm_win)==2 sm_win(1) specifies the height of the kernel
% in hertz and sm_win(2) the width of the kernel in seconds. Otherwise
% sm_win specifies the actual smoothing kernel.
%
% OUTPUTS:
%           C           --  complex valued time-frequency coherency
%                           [N,M] matrix with N frequencies and M
%                           time-points
%           F           --  frequency vector
%           T           --  time vector
%
% Note that due to the use of non-identical smoothing windows the resulting
% time-frequency coherency is not bound to [0,1].
% If no outputs are specified, time-frequency coherency is plotted.
%
% EXAMPLE:
% >>fs = 200; spec_win = fs; nfft = fs*3; tstep = fs/5;
% >>x1 = sin(2*pi*20*(1:fs*10)/fs); x2 = sin(2*pi*40*(1:fs*10)/fs);
% >>x = [x1,x1,x2]+randn(1,fs*30)/20; y = [x1,x2,x2]+randn(1,fs*30)/20;
% >>sm_win1 = [2,1.5]; sm_win2 = [100 5];
% >>tfcohf3(x,y,nfft,spec_win,sm_win1,sm_win2,tstep,1/fs)
%
% Time-frequency coherency between two signals sampled at 200 Hz of 30s
% duration for which synchronization jumps from 20 to 40 Hz. Data is
% decomposed using a 1s window. Cross spectrum is smoothed over an
% time-frequency area of 2Hz by 1.5s and the auto spectra over an
% time-frequency area of 100Hz by 5s.
%
% Please cite the following paper when using this code:
% Mehrkanoon S, Breakspear M, Daffertshofer A, Boonstra TW (2013). Non-
% identical smoothing operators for estimating time-frequency interdependence
% in electrophysiological recordings. EURASIP Journal on Advances in Signal
% Processing 2013, 2013:73. doi:10.1186/1687-6180-2013-73
%
% This work is licensed under a Creative Commons Attribution 3.0 Unported License
%
% T.W. Boonstra and S. Mehrkanoon          8-October-2012
% Systems Neuroscience Group, UNSW, Australia.
%
% See also FFT CONV

% if length(spec_win)==1
%     wl = spec_win;
% else
%     wl = length(spec_win);
% end
%
% % Zero-padding of signal
% x_new = zeros(length(x)+wl,1);
% y_new = x_new;
% x_new(fix(wl/2):fix(wl/2)+length(x)-1) = x;
% y_new(fix(wl/2):fix(wl/2)+length(x)-1) = y;

X = x_new;
Y = y_new;

L = size(X, 2);
fs = 0.5/dt;

% compute cross and auto spectra
% X = my_CWT(x_new,fscales,fs,Psi);
% Y = my_CWT(y_new,fscales,fs,Psi);

%XY = conj(Y) .* X;
XY = bsxfun(@times, Y, X);
X = abs(X).^2;
Y = abs(Y).^2;


% smooth cross spectrum using sm_win1

if(isscalar(time_winrad1))
    time_winrad1 = repmat(time_winrad1, [1 length(fscales)]);
end
if(length(time_winrad1) ~= length(fscales))
    error('Not a Valid Window Size!');
end

if(isscalar(scale_winrad1))
    scale_winrad1 = repmat(scale_winrad1, [1 length(fscales)]);
end
if(length(scale_winrad1) ~= length(fscales))
    error('Not a Valid Window Size!');
end

if (~isempty(time_winrad2))
    if(isscalar(time_winrad2))
        time_winrad2 = repmat(time_winrad2, [1 length(fscales)]);
    end
    
    if(length(time_winrad2) ~= length(fscales))
        error('Not a Valid Window Size!');
    end
    
    if(isscalar(scale_winrad2))
        scale_winrad2 = repmat(scale_winrad2, [1 length(fscales)]);
    end
    
    if(length(scale_winrad2) ~= length(fscales))
        error('Not a Valid Window Size!');
    end
end

XY_conv = complex(zeros(size(X,1), size(XY, 2)));

for f = 1:size(X,1)
    window = ones(2*scale_winrad1(f)+1,1) * ones(1,2*time_winrad1(f)+1);
    %     window = normpdf(-2*scale_winrad1(f):2*scale_winrad1(f),0,max(scale_winrad1(f),1)) * ...
    %              normpdf(-2*time_winrad1(f):2*time_winrad1(f),0,max(time_winrad1(f),1));
    window = window/sum(window(:));
    tmp_convolved = conv2(XY(max(f-scale_winrad1(f),1):min(f+scale_winrad1(f),end),:),window,'same');
    
    if(f <= size(X,1) - scale_winrad1(f))
        XY_conv(f,:) = tmp_convolved(min(f,scale_winrad1(f)+1),:);
    else
        XY_conv(f,:) = tmp_convolved(size(X,1) - f + 1,:);
    end
end

X_mean_power = zeros(size(X));
Y_mean_power = zeros(size(Y));

for f = 1:size(X,1)
    tmp_inds = find(coin(f,:)>0);
    if( isempty(time_winrad2) ||  (2*time_winrad2(f) + 1 > length(tmp_inds) - 1) )
        if(~isempty(tmp_inds))
            X_mean_power(f,:) = repmat(mean(X(f,max(tmp_inds(1) - time_winrad1(f),1):min(tmp_inds(end) + time_winrad1(f),end))),[1 size(X,2)]);
            Y_mean_power(f,:) = repmat(mean(Y(f,max(tmp_inds(1) - time_winrad1(f),1):min(tmp_inds(end) + time_winrad1(f),end))),[1 size(Y,2)]);
        end
        %         X_mean_power(f,:) = repmat(mean(X(f,:)),[1 size(X,2)]);
        %         Y_mean_power(f,:) = repmat(mean(Y(f,:)),[1 size(Y,2)]);
    else
        
        window = ones(2*scale_winrad2(f)+1,1) * ones(1,2*time_winrad2(f)+1);
        window = window/sum(window(:));
        
        tmp_convolved = conv2(X(max(f-scale_winrad2(f),1):min(f+scale_winrad2(f),end),:),window,'same');
        
        if(f <= size(X,1) - scale_winrad2(f))
            X_mean_power(f,:) = tmp_convolved(min(f,scale_winrad2(f)+1),:);
        else
            X_mean_power(f,:) = tmp_convolved(size(X,1) - f + 1,:);
        end
        
        tmp_convolved = conv2(Y(max(f-scale_winrad2(f),1):min(f+scale_winrad2(f),end),:),window,'same');
        
        if(f <= size(Y,1) - scale_winrad2(f))
            Y_mean_power(f,:) = tmp_convolved(min(f,scale_winrad2(f)+1),:);
        else
            Y_mean_power(f,:) = tmp_convolved(size(Y,1) - f + 1,:);
        end
        
        
    end
end


% compute tfcoh
Cxy = bsxfun(@rdivide, XY_conv, sqrt(bsxfun(@times, X_mean_power, Y_mean_power)));
%Cxy = XY_conv./sqrt(X_mean_power.*Y_mean_power);
%Cxy = Cxy .* coin;
Cxy = bsxfun(@times, Cxy, coin);

% if(~isempty(null_distribution))
%     for j=1:size(Cxy,1)
%         for i=1:size(Cxy,2)
%             Cxy(j,i) = length(find(Cxy(j,i) > null_distribution(:,j,i)))/size(null_distribution,1);
%         end
%     end
% end

% if no output arguments plot results
if nargout == 1
    varargout{1} = Cxy;
elseif nargout == 3;
    varargout{1} = Cxy;
    freq = scal2frq(fscales, wname, dt);%Fs*linspace(0,1,NFFT/2+1);
    varargout{2} = freq; %fs*linspace(0,1,length(select));%(select - 1)'*fs/nfft;
    varargout{3} = dt*[1:L];
else
    figure
    freq = scal2frq(fscales, wname, dt); %fs*linspace(0,1,length(select));%(select - 1)'*fs/nfft;
    time = dt*[1:length(x)];
    
    subplot(2,1,1)
    %wcoherence_plot(Cxy,freq);
    imagesc(time,freq,abs(Cxy))
    title('time-frequency coherency')
    xlabel('time [s]')
    ylabel('frequency [Hz]')
    
    subplot(2,1,2)
    imagesc(time,freq,tanh(abs(Cxy)))
    title('tanh of coherency')
    xlabel('time [s]')
    ylabel('frequency [Hz]')
end




function [coefs] = my_CWT(SIG,center_freqs, Fs, Psi)
dT = 0.5/Fs;
L = length(SIG);
L = L+mod(L,2);


coefs = zeros(length(center_freqs), length(SIG));

for curr_step=1:length(center_freqs)
    %     if (useWavToolbox)
    %         Psi = cmorwavf(-dT*L/2,dT*L/2,L,round(L/10),center_freqs(curr_step));
    %     else
    %Psi = cmplx_mor_wavelet(-dT*L/2,dT*L/2,L,round(L/10),center_freqs(curr_step));
    % end
    tmp = conv(Psi{curr_step}, SIG, 'same');
    %tmp = convFFT(Psi, SIG);
    
    %     if (curr_step == 1)
    %         coefs = zeros(length(center_freqs), length(tmp));
    %     end
    coefs(curr_step, :) = tmp;
end



function Psi = computeWavelet(SIG, Fs, center_freqs)


dT = 0.5/Fs;
L = length(SIG);
L = L+mod(L,2);

Psi = cell(1, length(center_freqs));

for curr_step = 1:length(center_freqs)
    
    Psi{curr_step} = cmplx_mor_wavelet(-dT*L/2,dT*L/2,L,round(L/10),center_freqs(curr_step));
    
end

function z = zscores(x)
%% Convert to z-scores

mu = mean(x, 2);
sigma = std(x, [], 2);
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus, x, mu);
z = bsxfun(@rdivide, z, sigma0);



function varargout = getcoi(scales, LSig, x)
%% Get cone of influence

bounds = [-8, 8];

% Compute the limits of cones of influence.
xR = LSig-bounds(2)*scales(end);
xL = -bounds(1)*scales(end)+1;
Lmin = [0 scales(end);xL 0];
Rmax = [xR 0;LSig scales(end)];

% Compute the cone of influence for each x value.
if ~isempty(x)
    x = x(:);
    scales = scales(:)';
    nbS = length(scales);
    nbX = length(x);
    L = repmat(bounds(1)*scales,nbX,1) + repmat(x,1,nbS);
    R = repmat(bounds(2)*scales,nbX,1) + repmat(x,1,nbS);
    
    PL = zeros(nbX,2);
    PR = zeros(nbX,2);
    for k = 1:nbX
        PL(k,:) = getLine([L(k,[1 end])' , scales([1 end])']);
        PR(k,:) = getLine([R(k,[1 end])' , scales([1 end])']);
    end
    
    cone = cell(1,nbX);
    warnState = warning;
    warning('Off'); %#ok<*WNOFF>
    for k = 1:nbX
        coneofx = zeros(nbS,LSig);
        L(k,:) = max(L(k,:),1);
        R(k,:) = min(R(k,:),LSig);
        for j=1:nbS
            coneofx(j,L(k,j):R(k,j)) = 1;
        end
        cone{k} = coneofx;
    end
    warning(warnState);
    if nbX==1 , cone = cone{1}; end
else
    L = []; R = []; cone = [];
end
Pmin = getLine(Lmin);
Pmax = getLine(Rmax);
xMin = (1-Pmin(2))/Pmin(1);
xMax = (1-Pmax(2))/Pmax(1);
if ~isempty(x)
    varargout = {cone,PL,PR,Pmin,Pmax,xMin,xMax};
else
    varargout = {Pmin,Pmax,xMin,xMax};
end

function P = getLine(L)

dx = L(2,1)-L(1,1);
dy = L(2,2)-L(1,2);
P(1) = dy/dx;
dxy = L(1,2)*L(2,1)-L(1,1)*L(2,2);
P(2) = dxy/dx;

function psi = cmplx_mor_wavelet(LB, UB, N, Fb, Fc)
%% Complex Morlet wavelet

X = linspace(LB,UB,N);
%psi = ((pi*Fb)^(-0.5))*exp(2*j*pi*Fc*X).*exp(-(X.*X)/Fb);
psi = ((pi*Fb)^(-0.5))*bsxfun(@times, exp(2*j*pi*Fc*X), exp(-bsxfun(@times, X, X)/Fb));


function z2 = convFFT(x,y)
n = length(x) + length(y) - 1;  % we need to zero-pad
z2 = ifft(bsxfun(@times, fft(x,n), fft(y,n)));
Nx = length(z2);
Nx = ceil(Nx/4);
z2 = z2(Nx + 1:Nx + length(x));
