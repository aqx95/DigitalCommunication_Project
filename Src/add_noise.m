% ------ Noise Addition ------
function noise_signal = add_noise(s0,signal_power,snr_db)
    n = randn(size(s0)); % generate noise samples
    
    SNR = power(10,snr_db/10); % calculate noise variance
    Np = signal_power/SNR; %calculate noise power
    
    noise = sqrt(Np/2).*n;
    noise_signal = s0 + noise;
end