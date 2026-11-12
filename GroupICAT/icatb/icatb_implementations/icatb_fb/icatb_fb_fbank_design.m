function [all_b,all_a,all_Ws,all_n] = icatb_fb_fbank_design(bandpass_ratio,fbank_opt)
%% design a filter bank arrays. This function uses matlab fucnctions to 
% estimate optimal filter order and filter Transfer function coefficients
%
% Inputs:
% 1. bandpass_ratio - a vector used to define filters cutoff frequencies
% 2. fbank_opt      - a struct used to define filters characteristics
%   fbank_opt
%	.Ws_value  : Filter transition width
%	.Rp        : Passband ripple expressed in dB (should be scalar)
%	.Rs        : Stopband ripple expressed in dB (should be scalar)
%	.time_n    : Number of time points (impacts only diagnostic plots)
%	.fft_n     : Number of fft points (impacts only diagnostic plots)
%	.filttype  : Filter type. choices are butter, cheby1, cheby2, and ellip
%	.Tr        : Sampling rate of the data
%	.plot_flag : If 1, the function will generate a diagnostic plot
%	.db_flag   : If 1, the plots will be in dB scale
%
% Outputs:
% all_b  - Transfer function numerator coefficients
% all_a  - Transfer function Denominator coefficients 
% all_Ws - Designed filters cutoff frequencies
% all_n  - Optimal filter order
%
% Example 1: running fbank_design without any inputs will return default
% values for fbank_opt
% fbank_opt = fbank_design;
%
% Example 2: designing 2 filters that divide the spectrum into 2 piece with
% the same length
% [all_b, all_a] = fbank_design([1,1],fbank_opt);
%%
switch nargin
    case 0
        %% return default fbank_opt
        fbank_opt.Ws_value = 0.05;
        fbank_opt.Rp = 3;
        fbank_opt.Rs = 30;
        fbank_opt.time_n = 30;
        fbank_opt.fft_n = 128;
        fbank_opt.filttype = 'cheby2';
        fbank_opt.Tr = 2;
        fbank_opt.plot_flag = 1;
        fbank_opt.db_flag = 1;
        all_b = fbank_opt;
        return 
    case 1
        %% read default fbank_opt and desing the fbanks
        fbank_opt.Ws_value = 0.05;
        fbank_opt.Rp = 3;
        fbank_opt.Rs = 30;
        fbank_opt.time_n = 30;
        fbank_opt.fft_n = 128;
        fbank_opt.filttype = 'cheby2';
        fbank_opt.Tr = 2;
        fbank_opt.plot_flag = 1;
        fbank_opt.db_flag = 1;
    case 2
    otherwise
        error('Wrong number of inputs.')
end

tmp_size = size(bandpass_ratio);
if tmp_size(1)==1
elseif tmp_size(2)==1
    bandpass_ratio = bandpass_ratio';
else
    error('bandpass_ratio should be a vector.')
end
    
Ws_value    = fbank_opt.Ws_value;
Rp          = fbank_opt.Rp;
Rs          = fbank_opt.Rs;
time_n      = fbank_opt.time_n;
fft_n       = fbank_opt.fft_n;
filttype    = fbank_opt.filttype;
Tr          = fbank_opt.Tr;
plot_flag   = fbank_opt.plot_flag;
db_flag     = fbank_opt.db_flag;

filter_n = length(bandpass_ratio);

bandpass_ratio  = bandpass_ratio/sum(bandpass_ratio);
bandpass_ratio = cumsum(bandpass_ratio);

for fn=1:filter_n
    if fn==1
        Wp = bandpass_ratio(1);
        Ws = bandpass_ratio(1) + Ws_value;
        
        if strcmp(filttype,'butter')
            [all_n{fn},Ws] = buttord(Wp,Ws,Rp,Rs);
            [b,a] = butter(all_n{fn},Ws);
        elseif strcmp(filttype,'cheby2')
            [all_n{fn},Ws] = cheb2ord(Wp,Ws,Rp,Rs);
            [b,a] = cheby2(all_n{fn},Rs,Ws);
        elseif strcmp(filttype,'cheby1')
            [all_n{fn},Ws] = cheb1ord(Wp,Ws,Rp,Rs);
            [b,a] = cheby1(all_n{fn},Rp,Wp);
        elseif strcmp(filttype,'ellip')
            [all_n{fn},Wp] = ellipord(Wp,Ws,Rp,Rs);
            [b,a] = ellip(all_n{fn},Rp,Rs,Wp);
        end      
        
        [all_hf{fn},w] = freqz(b,a,fft_n,'whole');
        [all_ht{fn},t] = impz(b,a,time_n);

        all_b{fn} = b;
        all_a{fn} = a;
        all_Ws(fn,:) = [0 Ws.*(1/Tr/2)];
        clear Wp Ws a b
    elseif fn==filter_n
        Wp = bandpass_ratio(fn-1);
        Ws = Wp + -1 * Ws_value;
        
        if strcmp(filttype,'butter')
            [all_n{fn},Ws] = buttord(Wp,Ws,Rp,Rs);
            [b,a] = butter(all_n{fn},Ws,'high');
        elseif strcmp(filttype,'cheby2')
            [all_n{fn},Ws] = cheb2ord(Wp,Ws,Rp,Rs);
            [b,a] = cheby2(all_n{fn},Rs,Ws,'high');
        elseif strcmp(filttype,'cheby1')
            [all_n{fn},Ws] = cheb1ord(Wp,Ws,Rp,Rs);
            [b,a] = cheby1(all_n{fn},Rp,Wp,'high');
        elseif strcmp(filttype,'ellip')
            [all_n{fn},Wp] = ellipord(Wp,Ws,Rp,Rs);
            [b,a] = ellip(all_n{fn},Rp,Rs,Wp,'high');
        end
        
%         [all_n{fn},Ws] = cheb2ord(Wp,Ws,Rp,Rs);
%         [b,a] = cheby2(all_n{fn},Rs,Ws,'high');
%         [b,a] = cheby2(fixed_n,Rs,Ws,'high');

        [all_hf{fn},w] = freqz(b,a,fft_n,'whole');
        [all_ht{fn},t] = impz(b,a,time_n);

        all_b{fn} = b;
        all_a{fn} = a;
        all_Ws(fn,:) = [Ws.*(1/Tr/2) (1/Tr/2)];
        clear Wp Ws a b
    else
        Wp = bandpass_ratio(fn-1:fn);
        Ws = Wp + [-1 1]*Ws_value;
        
        if strcmp(filttype,'butter')
            [all_n{fn},Ws] = buttord(Wp,Ws,Rp,Rs);
            [b,a] = butter(all_n{fn},Ws);
        elseif strcmp(filttype,'cheby2')
            [all_n{fn},Ws] = cheb2ord(Wp,Ws,Rp,Rs);
            [b,a] = cheby2(all_n{fn},Rs,Ws);
        elseif strcmp(filttype,'cheby1')
            [all_n{fn},Ws] = cheb1ord(Wp,Ws,Rp,Rs);
            [b,a] = cheby1(all_n{fn},Rp,Wp);
        elseif strcmp(filttype,'ellip')
            [all_n{fn},Wp] = ellipord(Wp,Ws,Rp,Rs);
            [b,a] = ellip(all_n{fn},Rp,Rs,Wp,'bandpass');
        end
        
        
%         [all_n{fn},Ws] = cheb2ord(Wp,Ws,Rp,Rs);
%         [b,a] = cheby2(all_n{fn},Rs,Ws);
%         [b,a] = cheby2(fixed_n,Rs,Ws);

        [all_hf{fn},w] = freqz(b,a,fft_n,'whole');
        [all_ht{fn},t] = impz(b,a,time_n);

        all_b{fn} = b;
        all_a{fn} = a;
        all_Ws(fn,:) = [Ws.*(1/Tr/2)];
        clear Wp Ws a b
    end
end

if plot_flag
    %%
    all_f = 1/Tr/2*linspace(0,1,fft_n/2+1);
    
%     tmp_fft_idx_tick = 1:3:median(1:fft_n)+1;
    tmp_fft_idx_tick = linspace(1,median(1:fft_n)+1,10);
    
    fft_idx_ticks = tmp_fft_idx_tick .* (1/Tr/fft_n);
    
    fft_idx_ticks = sort(all_Ws(:));
    idx_ticks_diff = [diff(fft_idx_ticks)']./max(fft_idx_ticks);
    idx_ticks_diff = idx_ticks_diff<1e-2;
    idx_ticks_diff = [0 idx_ticks_diff];
    fft_idx_ticks(find(idx_ticks_diff)) = [];
    
    fft_idx_tick_label = round((fft_idx_ticks*1000))/1000;
    
    
    filt_abs = abs(cell2mat(all_hf));
    ylab_txt  = 'Magnitude';
    
    fft_abs_min = min(filt_abs(:)).*1.1;
    fft_abs_min = max(fft_abs_min,-100);
    fft_abs_max = 1.2;
    
    if db_flag
        filt_abs = 20 .* log10(filt_abs);
        ylab_txt = 'Magnitude (dB)';
        
        fft_abs_min = min(filt_abs(:)).*1.1;
        fft_abs_min = max(fft_abs_min,-80);
        fft_abs_max = 1;
    end
    
    subplot(311)
    plot(all_f,filt_abs(1:fft_n/2+1,:),'LineWidth',2)
%     plot(all_f,log10(filt_abs(1:fft_n/2+1,:)) ,'LineWidth',2)

    ylim([fft_abs_min  fft_abs_max])
    xticks(fft_idx_ticks)
    xticklabels(fft_idx_tick_label)
    xtickangle(45)
    xlabel('Freq (Hz)')
    ylabel(ylab_txt)
    title('Frequency response of digital filter')
    
    subplot(312)
    fft_angle = angle(cell2mat(all_hf));
    fft_angle = fft_angle .* (abs(cell2mat(all_hf))>.5);
    plot(all_f,fft_angle(1:fft_n/2+1,:),'LineWidth',2)
    xticks(fft_idx_ticks)
    xticklabels(fft_idx_tick_label)
    xtickangle(45)
    xlabel('Freq (Hz)')
    ylabel('<H')
    title('Frequency response of digital filter')
    
    subplot(313)
    all_ht_mat = cell2mat(all_ht);
    all_ht_mat = all_ht_mat./max(all_ht_mat(:));
    hold on
    for fn=1:filter_n
        tmp_time_idx = (1:time_n) + floor((fn-1e-10)/4)*(time_n+round(time_n/5));
        tmp_bias = rem(fn-1,4);
        plot(tmp_time_idx,zeros(length(tmp_time_idx),1)+tmp_bias,'k','LineWidth',1)
        plot(tmp_time_idx,all_ht_mat(:,fn)+tmp_bias,'LineWidth',2)
    end
    hold off
    xticks([0 time_n])
    xlabel('Time')
    ylim([-1.5 4.5])
    title('Impulse response of digital filter')
end
end
