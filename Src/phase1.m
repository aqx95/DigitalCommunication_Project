%% Bits generator

close all
clear all

N = 1024; % no. of bits
snr = 0:5:50; % snr 0 to 50 dB in multiples of 5
signal_power = 1; 

% ----- Generate Binary Data -----
bin_data = binary_data(N); 


% ----- Transmission Data -----
trans_data = transmission_data(bin_data); 


% ----- Compute Error -----
error_list = [];
steps = 10; 
for k = 1: size(snr,2)    % traverse across col
    avg_error = 0;   %get average error
    for m = 1:steps
        noise_signal = add_noise(trans_data, signal_power, snr(k));
        received_signal = signal_threshold(noise_signal, 0);
        error_rate = compute_error(bin_data,received_signal);
        avg_error = avg_error + error_rate;
    end
    avg_error = avg_error / steps;
    error_list = [error_list,avg_error];
end

% ------------- Plot SNR vs. BER ----------------
plot(snr,error_list);
xlabel('Signal-to-Noise Ratio (SNR in dB)');
ylabel('Bit Error Rate');
title('SNR vs. BER');

semilogy (snr, error_list,'b*');
ylabel('Pe');
xlabel('Eb/No');
xlim([0 50]);
%-----------------------------------------------------------
