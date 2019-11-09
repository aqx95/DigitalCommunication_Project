
N = 1024;
br = 1000;       % bit rate of 1kps
bp= 1/br;            % bit period
L=16; %oversampling factor,L=Tb/Ts(Tb=bit period,Ts=sampling period)
%if a carrier is used, use L = Fs/Fc, where Fs >> 2xFc

Fc=10000; %carrier frequency
Fs=L*Fc;%sampling frequency

signal_length = Fs/br; 

bpsk_orig = binary_data(N);
%bpsk_orig = bpsk_binary_data;
bpsk_orig = 2.*bpsk_orig - 1;

% Get bits with signal length

bpsk_orig_signal = replicate(bpsk_orig, signal_length);

t1=bp/signal_length:bp/signal_length:length(bpsk_orig)*(bp);

subplot(5,1,1);
plot(t1,bpsk_orig_signal,'lineWidth',1.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Data Waveform');

% Carrier Wave
carrier = cos(2*pi*Fc*t1);
% subplot(5,1,2);
% plot(t1,carrier);grid on;
% axis([ 0 bp*10 -1.5 1.5]);
% ylabel('Amplitude(volt)');
% xlabel('Time(sec)');
% title('Carrier Signal');

% Modulation
bpsk_mod_signal = bpsk_orig_signal.*carrier;
subplot(5,1,2);
plot(t1,bpsk_mod_signal,'lineWidth',1.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Modulated Signal');

% Add Noise
bpsk_mod_noise_signal = add_noise(bpsk_mod_signal,1,10);
subplot(5,1,3);
plot(t1,bpsk_mod_noise_signal(1,:),'lineWidth',1.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Received Signal (Modulated Signal with Noise)');

% Demodulation
bpsk_demod_signal = bpsk_mod_noise_signal(1,:).*carrier;
subplot(5,1,4);
plot(t1,bpsk_demod_signal(1,:),'lineWidth',1.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Demodulated Signal');

int_op = [];
for i = 0:signal_length:length(bpsk_demod_signal)-signal_length
    int_ = trapz(bpsk_demod_signal(i+(signal_length/2):i+(signal_length/2)+1));
    int_op = [int_op int_];
end

%Threshold
th = 0;
bpsk_rec_bits = (round(int_op,1)> th);

bpsk_error_rate = compute_error(bpsk_orig, bpsk_rec_bits);

bpsk_rec_signal = replicate(bpsk_rec_bits, signal_length);
subplot(5,1,5);
plot(t1,bpsk_rec_signal,'lineWidth',1.5);grid on;
axis([ 0 bp*10 -.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Decoded Signal');

