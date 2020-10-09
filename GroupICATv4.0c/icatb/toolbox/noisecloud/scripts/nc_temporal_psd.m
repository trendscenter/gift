function feature = nc_temporal_psd(T)
% Calculate FFT of timeseries;
F = abs(fft(T,length(T)));
% Feature 24 to end: Power Spectrum Density, all frequencies;
try
    Hpsd = dspdata.psd(F(1:83),'Spectrumtype','Onesided');
    feature = Hpsd.Data';
catch
    feature = F(1:83);
end