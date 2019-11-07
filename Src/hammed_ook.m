% ------ OOK with error control codes -------

function ook_error_rate = ook(snr, ook_binary_data)
    
N = 1024;
br = 1000;       % bit rate of 1kps
bp= 1/br;            % bit period
L=16; %oversampling factor,L=Tb/Ts(Tb=bit period,Ts=sampling period)
%if a carrier is used, use L = Fs/Fc, where Fs >> 2xFc

Fc=10000; %carrier frequency
Fs=L*Fc; %sampling frequency

signal_length = Fs/br; 


%ook_orig = binary_data(N);
ook_orig = ook_binary_data;

ook_ham = encode(ook_orig,7,4,'hamming/fmt');

% Get bits with signal length
ook_orig_signal = replicate(ook_ham, signal_length);

t1=bp/signal_length:bp/signal_length:length(ook_ham)*(bp);
%subplot(5,1,1);
%plot(t1,ook_orig_signal,'lineWidth',2.5);grid on;
%axis([ 0 bp*10 -.5 1.5]);
%ylabel('Amplitude(volt)');
%xlabel('Time(sec)%');
%title('Transmitting signal');

% Carrier Wave
carrier = cos(2*pi*Fc*t1);
%subplot(5,1,2);
%plot(t1,carrier,'lineWidth',2.5);grid on;
%axis([ 0 bp*10 -1.5 1.5]);
%ylabel('Amplitude(volt)');
%xlabel(' Time(sec)');
%title('Carrier Signal');

% Modulation
ook_mod_signal = ook_orig_signal.*carrier;
%subplot(5,1,3);
%plot(t1,ook_mod_signal,'lineWidth',2.5);grid on;
%axis([ 0 bp*10 -1.5 1.5]);
%ylabel('Amplitude(volt)');
%xlabel(' Time(sec)');
%title('Modulated Signal');

% Add Noise
ook_mod_noise_signal = add_noise(ook_mod_signal,1,snr);
%subplot(5,1,4);
%plot(t1,ook_mod_noise_signal(1,:),'lineWidth',2.5);grid on;
%axis([ 0 bp*10 -1.5 1.5]);
%ylabel('Amplitude(volt)');
%xlabel('Time(sec)');
%title('Modulated Signal with Noise');

% Demodulation
ook_demod_signal = ook_mod_noise_signal(1,:).*carrier;
int_op = [];
for i = 0:signal_length:length(ook_demod_signal)-signal_length
    int_ = trapz(ook_demod_signal(i+(signal_length/2):i+(signal_length/2)+1));
    int_op = [int_op int_];
end
 
% Threshold
th = 0.5;
ook_rec_bits = (round(int_op,1)> th);
ook_dec = decode(ook_rec_bits,7,4,'hamming/fmt');

ook_rec_signal = replicate(ook_dec, signal_length);
%subplot(5,1,5);
%plot(t1,ook_rec_signal,'lineWidth',2.5);grid on;
%axis([ 0 bp*10 -.5 1.5]);
%ylabel('Amplitude(volt)');
%xlabel(' Time(sec)');
%title('Received Signal');

ook_error_rate = compute_error(ook_orig, ook_dec);

end