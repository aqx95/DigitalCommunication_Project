
% BFSK

N = 1024; 
br = 1000;      %bit rate
bp= 1/br;       %bit period                                             
L = 16; %oversampling factor

Fc = 10000; %carrier frequency
Fs = 16*Fc;  %sampling frequency

signal_length = Fs/br;  %length of signal

bfsk_orig = binary_data(N); %get original bit data
%bfsk_orig = bfsk_binary_data;

% Get bits with signal length
bfsk_orig_signal = replicate(bfsk_orig, signal_length);

t1=bp/signal_length:bp/signal_length:length(bfsk_orig)*(bp);  %time scaling

subplot(6,1,1);
plot(t1,bfsk_orig_signal,'lineWidth',1.5);grid on;
axis([ 0 bp*10 -.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Data Waveform Signal');

% Carrie5 Wave
f1 = Fc;
f2 = Fc/2;
t2 = bp/signal_length:bp/signal_length:bp; 
carrier1 = cos(2*pi*f1*t2);
carrier2 = cos(2*pi*f2*t2);

% Modulation               
bfsk_mod_signal=[];
for (i=1:1:length(bfsk_orig))
    if (bfsk_orig(i)==1)
        y = carrier1;
    else
        y = carrier2;
    end
    bfsk_mod_signal=[bfsk_mod_signal y];
end
subplot(6,1,2);
plot(t1,bfsk_mod_signal,'lineWidth',1.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
title('Modulated Signal');

% Add Noise
bfsk_noise_signal = add_noise(bfsk_mod_signal,1,10); 
subplot(6,1,3);
plot(t1,bfsk_noise_signal,'lineWidth',1.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
title('Received Signal (Modulated Signal with Noise)');

% Demodulation
y1 = bfsk_noise_signal.*cos(2*pi*f1*t1);
y2 = bfsk_noise_signal.*cos(2*pi*f2*t1);

subplot(6,1,4);
plot(t1,y1,'lineWidth',1.5);grid on;
axis([ 0 bp*10 -1.5 1.5]); 
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
title('Demodulated Signal (f1)');

subplot(6,1,5);
plot(t1,y2,'lineWidth',1.5);grid on;
axis([ 0 bp*10 -1.5 1.5]); 
xlabel('Time(sec)');
ylabel('Amplitude(volt)');
title('Demodulated Signal (f2)');


int1_bfsk = [];
for i = 0:signal_length:length(y1)-signal_length
    int_ = trapz(y1(i+(signal_length/2):i+(signal_length/2)+1));
    int1_bfsk = [int1_bfsk int_];
end
int2_bfsk = [];
for i = 0:signal_length:length(y2)-signal_length
    int_ = trapz(y2(i+(signal_length/2):i+(signal_length/2)+1));
    int2_bfsk = [int2_bfsk int_];
end
int = int1_bfsk - int2_bfsk;

bfsk_demod_signal = signal_threshold(int, 0);

bfsk_error_rate = compute_error(bfsk_orig, bfsk_demod_signal);

bfsk_demod_signal = replicate(bfsk_demod_signal, signal_length);

subplot(6,1,6)
plot(t1,bfsk_demod_signal,'lineWidth',1.5);grid on;
axis([0 bp*10 -.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Decoded Signal')


