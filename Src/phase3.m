%% Plot error of OOK and OOK with ham
clear all
close all

N = 1024;

ook_binary_data = binary_data(N);

snr = 0:2:20;
ber = (1:11);
steps = 11; 
for k = 1:11
    avg_error = 0;   %get average error
    er = (1:10);
    for m = 1:steps
        er(m) = ook(snr(k), ook_binary_data);
        avg_error = avg_error + er(m);
    end
    avg_error = avg_error / steps;
    ber(k) = avg_error;
end

ber_ham = (1:11);
for k = 1:11
   avg_error = 0;   %get average error
   er = (1:10);
    for m = 1:steps
        er(m) = hammed_ook(snr(k), ook_binary_data);
        avg_error = avg_error + er(m);
    end
    avg_error = avg_error / steps;
    ber_ham(k) = avg_error;
end

semilogy(snr,ber,'*')
hold on
semilogy(snr,ber)
hold on
semilogy(snr,ber_ham,'*')
hold on
semilogy(snr,ber_ham)
grid
legend('OOK','','hammed OOK','')
xlabel('Eb/No (dB)')
ylabel('Bit Error Rate')
title('Bit Error Rate Performance of OOK')
xlim([0 18])