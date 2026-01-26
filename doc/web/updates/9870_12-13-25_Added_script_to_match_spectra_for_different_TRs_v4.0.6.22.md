# Added script to match spectra for different TRs v4.0.6.22 <BR>
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
    
