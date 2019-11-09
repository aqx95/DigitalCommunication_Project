% global bin_data;
% global N;
% global ber, 
% global ber_ook;
% global ber_bpsk;
% global ber_bfsk;
% global snr;
% global steps;

% Plot error of OOK
clc
clear all
close all

N = 1024;
bin_data = binary_data(N);
snr = 0:2:20;
ber = (1:11);
steps = 11; 

for k = 1:11
    avg_error = 0;   %get average error
    er = (1:10);
    for m = 1:steps
        er(m) = ook(snr(k), bin_data);
        avg_error = avg_error + er(m);
    end
    avg_error = avg_error / steps;
    ber(k) = avg_error;
end

% semilogy(snr,ber,'*')
% hold on
% semilogy(snr,ber)
% grid
% legend('OOK','')
% xlabel('Eb/No (dB)')
% ylabel('Bit Error Rate')
% title('Bit Error Rate Performance of OOK')
% xlim([0 15])

ber_ook = ber;

% Plot error of BFSK
% clc
% clear all
% close all

N = 1024;
%bin_data = binary_data(N);
snr = 0:2:20;
ber = (1:11);
steps = 11; 

for k = 1:11
    avg_error = 0;   %get average error
    er = (1:10);
    for m = 1:steps
        er(m) = bfsk(snr(k), bin_data);
        avg_error = avg_error + er(m);
    end
    avg_error = avg_error / steps;
    ber(k) = avg_error;
end

% semilogy(snr,ber,'*')
% hold on
% semilogy(snr,ber)
% grid
% legend('BFSK','')
% xlabel('Eb/No (dB)')
% ylabel('Bit Error Rate')
% title('Bit Error Rate Performance of BFSK')
% xlim([0 15])

ber_bfsk = ber;

% Plot error of BPSK
% clc
% clear all
% close all

N = 1024;
bin_data = binary_data(N);
snr = 0:2:20;
ber = (1:11);
steps = 11; 

for k = 1:11
    avg_error = 0;   %get average error
    er = (1:10);
    for m = 1:steps
        er(m) = bpsk(snr(k), bin_data);
        avg_error = avg_error + er(m);
    end
    avg_error = avg_error / steps;
    ber(k) = avg_error;
end

% semilogy(snr,ber,'*')
% hold on
% semilogy(snr,ber)
% grid
% legend('BPSK','')
% xlabel('Eb/No (dB)')
% ylabel('Bit Error Rate')
% title('Bit Error Rate Performance of BPSK')
% xlim([0 5])

ber_bpsk = ber;

% PLot all - run first 3 sections first before running this

% ook
semilogy(snr,ber_ook,'*')
hold on
semilogy(snr,ber_ook)
hold on 

% bpsk
semilogy(snr,ber_bpsk,'*')
hold on
semilogy(snr,ber_bpsk)
hold on 

% bfsk
semilogy(snr,ber_bfsk,'*')
hold on
semilogy(snr,ber_bfsk)
hold on 

grid
legend('OOK','','BPSK','','BFSK','')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
title('Bit Error Rate Performance of BPSK')
xlim([0 20])
