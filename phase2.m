%% BFSK


x=[1,0,1,1,0,0,1,0]                                   % Binary Information
bp= 0.001;                                                    % bit period
disp(' Binary information at Trans mitter :');
disp(x);
%XX representation of transmitting binary information as digital signal XXX
bit=[]; 
for n=1:1:length(x)
    if x(n)==1;
       se=ones(1,100);
    else x(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];
end
t1=bp/100:bp/100:length(x)*(bp);
subplot(3,1,1);
plot(t1,bit,'lineWidth',2.5);grid on;
axis([ 0 bp*length(x) -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('transmitting information as digital signal');
%XXXXXXXXXXXXXXXXXXXXXXX Binary-FSK modulation XXXXXXXXXXXXXXXXXXXXXXXXXXX%
A=1;                                          % Amplitude of carrier signal
br=1/bp;                                                         % bit rate
f1=br*8;                           % carrier frequency for information as 1
f2=br*2;                           % carrier frequency for information as 0
t2=bp/100:bp/100:bp;                 
ss=length(t2);
m=[];
for (i=1:1:length(x))
    if (x(i)==1)
        y=A*cos(2*pi*f1*t2);
    else
        y=A*cos(2*pi*f2*t2);
    end
    m=[m y];
end
t3=bp/100:bp/100:100*length(x)*(bp/100);
subplot(3,1,2);
plot(t3,m);
xlabel('time(sec)');
ylabel('amplitude(volt)');
title('waveform for binary FSK modulation coresponding binary information');
%XXXXXXXXXXXXXXXXXXXX Binary FSK demodulation XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
mn=[];
for n=ss:ss:length(m)
  t=bp/100:bp/100:bp;
  y1=cos(2*pi*f1*t);                    % carrier siignal for information 1
  y2=cos(2*pi*f2*t);                    % carrier siignal for information 0
  mm=y1.*m((n-(ss-1)):n);
  mmm=y2.*m((n-(ss-1)):n);
  t4=bp/99:bp/99:bp;
  z1=trapz(t4,mm)                                             % intregation 
  z2=trapz(t4,mmm)                                            % intregation 
  zz1=round(2*z1/bp)
  zz2= round(2*z2/bp)
  if(zz1>A/2)      % logic lavel= (0+A)/2 or (A+0)/2 or 2.5 ( in this case)
    a=1;
  else(zz2>A/2)
    a=0;
  end
  mn=[mn a];
end
disp(' Binary information at Reciver :');
disp(mn);
%XXXXX Representation of binary information as digital signal which achived 
%after demodulation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
bit=[];
for n=1:length(mn);
    if mn(n)==1;
       se=ones(1,100);
    else mn(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];
end
t4=bp/100:bp/100:100*length(mn)*(bp/100);
subplot(3,1,3)
plot(t4,bit,'LineWidth',2.5);grid on;
axis([ 0 bp*length(mn) -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('recived information as digital signal after binary FSK demodulation');

%% OOK

br = 1000;       % bit rate of 1kps
bp= 1/br;            % bit period
L=16; %oversampling factor,L=Tb/Ts(Tb=bit period,Ts=sampling period)
%if a carrier is used, use L = Fs/Fc, where Fs >> 2xFc

Fc=10000; %carrier frequency
Fs=L*Fc;%sampling frequency


ook_orig = get_baseband_data(N);

% Duplicate bits to sync with bit period
ook_orig_signal = sync(ook_orig);

t1=bp/100:bp/100:length(ook_orig)*(bp);
subplot(5,1,1);
plot(t1,ook_orig_signal,'lineWidth',2.5);grid on;
axis([ 0 bp*10 -.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Transmitting signal');

% Carrier Wave
carrier = cos(2*pi*Fc*t1);
subplot(5,1,2);
plot(t1,carrier,'lineWidth',2.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Carrier Signal');

% Modulation
ook_mod_signal = ook_orig_signal.*carrier;
subplot(5,1,3);
plot(t1,ook_mod_signal,'lineWidth',2.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Modulated Signal');

% Add Noise
ook_mod_noise_signal = get_unipolar_received_data(ook_mod_signal);
subplot(5,1,4);
plot(t1,ook_mod_noise_signal(1,:),'lineWidth',2.5);grid on;
axis([ 0 bp*10 -1.5 1.5]);
ylabel('Amplitude(volt)');
xlabel('Time(sec)');
title('Modulated Signal with Noise');

% Demodulation
ook_demod_signal = ook_mod_noise_signal(1,:).*carrier;
int_op = [];
s = 100;
for i = 0:s:length(ook_demod_signal)-s
    int_ = trapz(ook_demod_signal(i+(s/2):i+(s/2)+1));
    int_op = [int_op int_];
end
 
% Threshold
th = 0.5;
ook_rec_bits = (round(int_op,1)>= th);
ook_rec_signal = sync(ook_rec_bits);
subplot(5,1,5);
plot(t1,ook_rec_signal,'lineWidth',2.5);grid on;
axis([ 0 bp*10 -.5 1.5]);
ylabel('Amplitude(volt)');
xlabel(' Time(sec)');
title('Received Signal');



%% BFSK


fc = 1e5;       % Carrier freq
fs = 16*fc;         % Sampling freq


% --------------- Duplicate bit info to syn with bit period --------------
bit=[]; 
for n=1:1:length(y0)
    if y0(n)==1;
       se=ones(1,100);
    else y0(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];
end
t1=bp/100:bp/100:length(y0)*(bp);
subplot(3,1,1);
plot(t1,bit,'lineWidth',2.5);grid on;
axis([ 0 bp*10 -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('transmitting information as digital signal');

% ---------------------- Binary-FSK modulation -----------------------------
A=1;                                          % Amplitude of carrier signal                                                        % bit rate
f1=br*8;                           % carrier frequency for information as 1
f2=br*2;                           % carrier frequency for information as 0
t2=bp/100:bp/100:bp;                 
ss=length(t2);
m=[];
for (i=1:1:length(y0))
    if (y0(i)==1)
        y=A*cos(2*pi*f1*t2);
    else
        y=A*cos(2*pi*f2*t2);
    end
    m=[m y];
end
t3=bp/100:bp/100:100*length(y0)*(bp/100);
subplot(3,1,2);
plot(t3,m);
axis([ 0 bp*10 -1.5 1.5]);
xlabel('Time(secs)');
ylabel('Amplitude(volts)');
title('Modulated Signal');

%------------------- Adding Additive WHite Noise ---------------------
n = randn(size(m)); % generate noise samples
nvar = power(10,-10/10); % calculate noise variance
m = m + nvar.*n;

%-------------------- Binary FSK demodulation ------------------------
mn=[];
for n=ss:ss:length(m)
  t=bp/100:bp/100:bp;
  y1=cos(2*pi*f1*t);                    % carrier siignal for information 1
  y2=cos(2*pi*f2*t);                    % carrier siignal for information 0
  mm=y1.*m((n-(ss-1)):n);
  mmm=y2.*m((n-(ss-1)):n);
  t4=1:1:length(mm);
  z1=trapz(t4,mm);                                             % intregation 
  z2=trapz(t4,mmm);                                            % intregation 
  zz1=round(2*z1/bp);
  zz2= round(2*z2/bp);
  if(zz1>A/2)      % logic lavel= (0+A)/2 or (A+0)/2 or 2.5 ( in this case)
    a=1;
  else(zz2>A/2)
    a=0;
  end
  mn=[mn a];
end
% disp(' Binary information at Reciver :');
% disp(mn);
%XXXXX Representation of binary information as digital signal which achived 
%after demodulation XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
bit=[];
for n=1:length(mn);
    if mn(n)==1;
       se=ones(1,100);
    else mn(n)==0;
        se=zeros(1,100);
    end
     bit=[bit se];
end
t4=bp/100:bp/100:100*length(mn)*(bp/100);
subplot(3,1,3)
plot(t4,bit,'LineWidth',2.5);grid on;
axis([ 0 bp*10 -.5 1.5]);
ylabel('amplitude(volt)');
xlabel(' time(sec)');
title('recived information as digital signal after binary FSK demodulation');




%% FUNCTIONS




% ----- returns error rate by comparing original and received signal ----
function error_rate = compute_error_rate(N,snr)
    % compute error rate
    error_rate = 0;
    steps = 20;
    snR = snr;
    avg_error = 0;
    for i = 1:steps
        Error = 0;
        orig = get_baseband_data(N);
        trans = 2.*orig - 1;
        noise_signal = add_noise(trans, snR);
        for j = 1:N
            if noise_signal(j)<0
                noise_signal(j) = -1;
            else
                noise_signal(j) = 1;
            end
        end
        for k = 1:1:N
            if trans(k) ~= noise_signal(k)
                Error = Error + 1;
            end
        end
        error_rate = Error/N;
        avg_error = avg_error + error_rate;
    end
    error_rate = avg_error/steps;
end

function sync_signal = sync(data_bits)
    ook_orig_signal = [];
    for n=1:1:length(data_bits)
        if data_bits(n)==1;
        se=ones(1,100);
        else data_bits(n)==0;
            se=zeros(1,100);
        end
        ook_orig_signal=[ook_orig_signal se];
    end
    sync_signal = ook_orig_signal;
end



