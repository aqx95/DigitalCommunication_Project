function bpsk_error_rate = bpsk(snr, bpsk_binary_data)

N = 1024;
br = 1000;       % bit rate of 1kps
bp= 1/br;            % bit period
L=16; %oversampling factor,L=Tb/Ts(Tb=bit period,Ts=sampling period)
%if a carrier is used, use L = Fs/Fc, where Fs >> 2xFc

Fc=10000; %carrier frequency
Fs=L*Fc;%sampling frequency

signal_length = Fs/br; 

bpsk_orig_binary = bpsk_binary_data;
bpsk_orig = zeros(1,N);
for m =1:length(bpsk_orig_binary)
    if bpsk_orig_binary(m) == 1
        bpsk_orig(m) = 1;
    else
        bpsk_orig(m) = -1;
    end
end

% Get bits with signal length
bpsk_orig_signal = replicate(bpsk_orig, signal_length);

t1=bp/signal_length:bp/signal_length:length(bpsk_orig)*(bp);

% subplot(6,1,1);
% plot(t1,bpsk_orig_signal,'lineWidth',2.5);grid on;
% axis([ 0 bp*10 -1.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel('Time(sec)');
% title('Transmitting signal');

% Carrier Wave
carrier = cos(2*pi*Fc*t1);
% subplot(6,1,2);
% plot(t1,carrier);grid on;
% axis([ 0 bp*10 -1.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Carrier Signal');

% Modulation
bpsk_mod_signal = bpsk_orig_signal.*carrier;
% subplot(6,1,3);
% plot(t1,bpsk_mod_signal);grid on;
% axis([ 0 bp*10 -1.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Modulated Signal');

% Add Noise
bpsk_mod_noise_signal = add_noise(bpsk_mod_signal,1,snr);
% subplot(6,1,4);
% plot(t1,bpsk_mod_noise_signal,'lineWidth',2.5);grid on;
% axis([ 0 bp*10 -1.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel('Time(sec)');
% title('Modulated Signal with Noise');

% Demodulation
bpsk_int_signal = bpsk_mod_noise_signal.*(2*carrier);
[b,a] = butter(6,0.2);
bpsk_demod_signal = filtfilt(b,a,bpsk_int_signal);
% subplot(6,1,5);
% plot(t1,bpsk_demod_signal,'lineWidth',2.5);grid on;
% axis([ 0 bp*10 -.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Demodulated Signal');

%Threshold
th = 0;
bpsk_rec_signal = zeros(1,N);
index = 1;
for i = signal_length/2:signal_length:length(bpsk_demod_signal)
    if bpsk_demod_signal(i) > th
        bpsk_rec_signal(index) = 1;
    else
        bpsk_rec_signal(index) = 0;
    end
    index = index + 1;
end
    
bpsk_error_rate = compute_error(bpsk_orig_binary, bpsk_rec_signal);
bpsk_rec_signal2 = replicate(bpsk_rec_signal, signal_length);

% subplot(6,1,6);
% plot(t1,bpsk_rec_signal2,'lineWidth',2.5);grid on;
% axis([ 0 bp*10 -.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel(' Time(sec)');
% title('Received Signal');

end
