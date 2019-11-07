N = 1024;
br = 1000;       % bit rate of 1kps
bp= 1/br;            % bit period
L=16; %oversampling factor,L=Tb/Ts(Tb=bit period,Ts=sampling period)
%if a carrier is used, use L = Fs/Fc, where Fs >> 2xFc

Fc=10000; %carrier frequency
Fs=L*Fc;%sampling frequency

signal_length = Fs/br; 

bpsk_orig = binary_data(N);
for k =1:length(bpsk_orig)
    if bpsk_orig(k) == 1
        bpsk_orig(k) = 1;
    else
    bpsk_orig(k)=-1;
    end
end

% Get bits with signal length

bpsk_orig_signal = replicate(bpsk_orig, signal_length);

t1=bp/signal_length:bp/signal_length:length(bpsk_orig)*(bp);

subplot(5,1,1);
plot(t1,bpsk_orig_signal,'lineWidth',2.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Transmitting signal');

% Carrier Wave
carrier = cos(2*pi*Fc*t1);
subplot(5,1,2);
plot(t1,carrier);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Carrier Signal');

% Modulation
bpsk_mod_signal = bpsk_orig_signal.*carrier;
subplot(5,1,3);
plot(t1,bpsk_mod_signal);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Modulated Signal');

% Add Noise
bpsk_mod_noise_signal = add_noise(bpsk_mod_signal,1,5);
subplot(5,1,4);
plot(t1,bpsk_mod_noise_signal(1,:),'lineWidth',2.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Modulated Signal with Noise');

% Demodulation
bpsk_demod_signal = bpsk_mod_signal(1,:).*carrier;
int_op = [];
for i = 0:signal_length:length(bpsk_demod_signal)-signal_length
    int_ = trapz(bpsk_demod_signal(i+(signal_length/2):i+(signal_length/2)+1));
    int_op = [int_op int_];
end

%Threshold
th = 0.5;
ook_rec_bits = (round(int_op,1)> th);
ook_rec_signal = replicate(ook_rec_bits, signal_length);
subplot(5,1,5);
plot(t1,ook_rec_signal,'lineWidth',2.5);grid on;
axis([ 0 bp*10 -.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Received Signal');


