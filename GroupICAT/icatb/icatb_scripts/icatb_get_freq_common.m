function [SA_common, SB_common] = icatb_get_freq_common(s_c1_file_A, s_c1_file_B, paramsA, paramsB)
    % Cyrus Eierud 12/13/2025 TReNDS
    % Function that creates equal sample space for data
    %  with different TRs
    %
    % Inputs:
    % s_c1_file_A = path and file to {prefix}_ica_c1-1.mat for first dataset
    % s_c1_file_B = path and file to {prefix}_ica_c1-1.mat for second dataset (with different TR)
    % paramsA.Fs  = Frequency of Imaging (1/TR)
    % paramsA.tapers = [3 5]; See icatb_mtspectrumc.m for more info
    % paramsA.pad    = 0; %Needs to be zero for icatb_get_freq_common (more info in icatb_mtspectrumc.m)
    % paramsA.fpass  = [0 1/(2*TR)]; % [fmin fmax] this fun supports fmin=0
    % paramsB is analogous with paramsA structure
    %
    % Outputs:
    % SA_common is the array containing the frequencies matching SB_common
    % SB_common is the array containing the frequencies matching SA_common
    % 
    % Example:
    % % Case A: TR=0.6667s
    % paramsA.Fs     = 1.5;
    % paramsA.tapers = [3 5];
    % paramsA.pad    = 0;
    % paramsA.fpass  = [0 0.75];
    % s_c1_file_A='/path/to/prefix_ica_c1-1.mat';
    % 
    % % Case B: TR=2s
    % paramsB.Fs     = 0.5;
    % paramsB.tapers = [3 5];
    % paramsB.pad    = 0;
    % paramsB.fpass  = [0 0.25];
    % s_c1_file_B='/different_path/to/prefix_different_ica_c1-1.mat';
    %
    % % Running the function to get spectra in matching freq below
    % [SA_common, SB_common] = icatb_get_freq_common(s_c1_file_A, s_c1_file_B, paramsA, paramsB)
    
    stru_A = load(s_c1_file_A);
    stru_B = load(s_c1_file_B);

    [SA, fA] = icatb_mtspectrumc(stru_A.tc, paramsA);
    [SB, fB] = icatb_mtspectrumc(stru_B.tc, paramsB);

    % Common frequency grid over the shared band from zero to the minimal max frequency (Hz)
    d_min_max_f = min([paramsA.fpass(2) paramsB.fpass(2)]);
    df_common = 0.001;            % choose your preferred spacing
    f_common  = (0:df_common:d_min_max_f)';

    % Interpolate (linear; use 'pchip' if you want smoother)
    SA_common = interp1(fA, SA, f_common, 'linear', 'extrap');
    SB_common = interp1(fB, SB, f_common, 'linear', 'extrap');

    figure();
    plot(f_common, SA_common(:,1));
    hold on;
    plot(f_common, SB_common(:,1));
    xlabel('Frequency [Hz]')
    ylabel('Power [arbitrary units]')
    title('Plotting first component for demonstration')
    legend('Spectra from Data A','Spectra from Data B');

end