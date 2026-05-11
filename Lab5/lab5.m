clear; clc; close all;

%% PARAMETERS

L1 = 2^13;
L2 = 2^12;

Q = 64;

fs = 30e6;
fc = 350e6;
M = 100;

M1 = 4;
M2 = 8;

rolloff = 1/3;
S = 10;

omega1 = -pi/3;
omega2 = pi/6;

wc = 2*pi*fc/(M*fs);

SNR_test = 5:20;
Neq_test = 0:10;

c_ideal = 1;
c = 0.25*exp(1i*0.1*pi)*[1 zeros(1,15) 2.4 zeros(1,15) 1];

%% SIGNALS

d1 = randi([0 Q-1],L1,1);
d2 = randi([0 Q-1],L2,1);

x1 = qammod(d1,Q,'UnitAveragePower',true).';
x2 = qammod(d2,Q,'UnitAveragePower',true).';

%% FILTERS

G1 = rcosdesign(rolloff,S,M1,'sqrt')/sqrt(M1);
G2 = rcosdesign(rolloff,S,M2,'sqrt')/sqrt(M2);

H1 = M1*G1;
H2 = M2*G2;

G_lp = rcosdesign(rolloff,S,M,'sqrt')/sqrt(M);
H_lp = M*G_lp;

%% TRANSMUX TRANSMITTER

x1_tx = upsample(x1,M1);
x1_tx = conv(x1_tx,H1);
n = 0:length(x1_tx)-1;
x1_tx = x1_tx .* exp(1i*omega1*n);

x2_tx = upsample(x2,M2);
x2_tx = conv(x2_tx,H2);
n = 0:length(x2_tx)-1;
x2_tx = x2_tx .* exp(1i*omega2*n);

L = max(length(x1_tx),length(x2_tx));

x1_tx = [x1_tx zeros(1,L-length(x1_tx))];
x2_tx = [x2_tx zeros(1,L-length(x2_tx))];

y_tx = x1_tx + x2_tx;

%% LAB 3 REFERENCE CASE

y_rx_lab3 = y_tx;

n = 0:length(y_rx_lab3)-1;

x1_rx = y_rx_lab3 .* exp(-1i*omega1*n);
x1_rx = conv(x1_rx,G1);
x1_rx = downsample(x1_rx,M1);
x1_lab3 = x1_rx(S+1:S+L1);

x2_rx = y_rx_lab3 .* exp(-1i*omega2*n);
x2_rx = conv(x2_rx,G2);
x2_rx = downsample(x2_rx,M2);
x2_lab3 = x2_rx(S+1:S+L2);

SIDR1_lab3 = 10*log10(sum(abs(x1).^2)/sum(abs(x1_lab3-x1).^2));
SIDR2_lab3 = 10*log10(sum(abs(x2).^2)/sum(abs(x2_lab3-x2).^2));

%% IDEAL ZERO-IF CASE

y = upsample(y_tx,M);
y = conv(y,H_lp);

m = 0:length(y)-1;
y = real(y .* exp(1i*wc*m));

y = conv(y,c_ideal);

m = 0:length(y)-1;
y = 2*y .* exp(-1i*wc*m);

y = conv(y,G_lp);
y = downsample(y,M);

y_rx_ideal = y(S+1:S+length(y_tx));

n = 0:length(y_rx_ideal)-1;

x1_rx = y_rx_ideal .* exp(-1i*omega1*n);
x1_rx = conv(x1_rx,G1);
x1_rx = downsample(x1_rx,M1);
x1_ideal = x1_rx(S+1:S+L1);

x2_rx = y_rx_ideal .* exp(-1i*omega2*n);
x2_rx = conv(x2_rx,G2);
x2_rx = downsample(x2_rx,M2);
x2_ideal = x2_rx(S+1:S+L2);

SIDR1_ideal = 10*log10(sum(abs(x1).^2)/sum(abs(x1_ideal-x1).^2));
SIDR2_ideal = 10*log10(sum(abs(x2).^2)/sum(abs(x2_ideal-x2).^2));

err1_ideal = sum(qamdemod(x1_ideal,Q,'UnitAveragePower',true).' ~= d1);
err2_ideal = sum(qamdemod(x2_ideal,Q,'UnitAveragePower',true).' ~= d2);

fprintf("Lab 3 reference case\n");
fprintf("SIDR1 = %.2f dB\n",SIDR1_lab3);
fprintf("SIDR2 = %.2f dB\n\n",SIDR2_lab3);

fprintf("Ideal zero-IF case\n");
fprintf("SIDR1 = %.2f dB, symbol errors x1 = %d\n",SIDR1_ideal,err1_ideal);
fprintf("SIDR2 = %.2f dB, symbol errors x2 = %d\n\n",SIDR2_ideal,err2_ideal);

%% SPECTRUM PLOT

NFFT = 2^nextpow2(length(y_tx));

YTX = fftshift(fft(y_tx,NFFT));
YRX = fftshift(fft(y_rx_ideal,NFFT));

f = linspace(-fs/2,fs/2,NFFT);

figure;
plot(f/1e6,20*log10(abs(YTX)/max(abs(YTX)) + eps)); hold on;
plot(f/1e6,20*log10(abs(YRX)/max(abs(YRX)) + eps));
grid on;
xlabel('Frequency [MHz]');
ylabel('Magnitude [dB]');
legend('y_{tx}','y_{rx} ideal');
title('Spectrum before and after ideal zero-IF system');

%% CHANNEL WITHOUT EQUALIZER

y = upsample(y_tx,M);
y = conv(y,H_lp);

m = 0:length(y)-1;
y = real(y .* exp(1i*wc*m));

y = conv(y,c);

m = 0:length(y)-1;
y = 2*y .* exp(-1i*wc*m);

y = conv(y,G_lp);
y = downsample(y,M);

y_rx_channel = y(S+1:S+length(y_tx));

n = 0:length(y_rx_channel)-1;

x1_rx = y_rx_channel .* exp(-1i*omega1*n);
x1_rx = conv(x1_rx,G1);
x1_rx = downsample(x1_rx,M1);
x1_channel = x1_rx(S+1:S+L1);

x2_rx = y_rx_channel .* exp(-1i*omega2*n);
x2_rx = conv(x2_rx,G2);
x2_rx = downsample(x2_rx,M2);
x2_channel = x2_rx(S+1:S+L2);

SINDR1_channel = 10*log10(sum(abs(x1).^2)/sum(abs(x1_channel-x1).^2));
SINDR2_channel = 10*log10(sum(abs(x2).^2)/sum(abs(x2_channel-x2).^2));

err1_channel = sum(qamdemod(x1_channel,Q,'UnitAveragePower',true).' ~= d1);
err2_channel = sum(qamdemod(x2_channel,Q,'UnitAveragePower',true).' ~= d2);

fprintf("Channel without equalizer\n");
fprintf("SINDR1 = %.2f dB, symbol errors x1 = %d\n",SINDR1_channel,err1_channel);
fprintf("SINDR2 = %.2f dB, symbol errors x2 = %d\n\n",SINDR2_channel,err2_channel);

%% EQUALIZATION FOR ALL SNR AND Neq VALUES

SINDR1 = zeros(length(Neq_test),length(SNR_test));
SINDR2 = zeros(length(Neq_test),length(SNR_test));

symerr1 = zeros(length(Neq_test),length(SNR_test));
symerr2 = zeros(length(Neq_test),length(SNR_test));

for snr_idx = 1:length(SNR_test)

    SNR = SNR_test(snr_idx);

    y = upsample(y_tx,M);
    y = conv(y,H_lp);

    m = 0:length(y)-1;
    y = real(y .* exp(1i*wc*m));

    y = conv(y,c);
    y = awgn(y,SNR,'measured');

    m = 0:length(y)-1;
    y = 2*y .* exp(-1i*wc*m);

    y = conv(y,G_lp);
    y = downsample(y,M);

    y_rx = y(S+1:end);

    for neq_idx = 1:length(Neq_test)

        Neq = Neq_test(neq_idx);

        [~,d_min,~,y_rx_eq] = fir_eq(y_tx,y_rx,Neq);

        y_rx_eq = y_rx_eq(d_min+1:d_min+length(y_tx));

        n = 0:length(y_rx_eq)-1;

        x1_rx = y_rx_eq .* exp(-1i*omega1*n);
        x1_rx = conv(x1_rx,G1);
        x1_rx = downsample(x1_rx,M1);
        x1_est = x1_rx(S+1:S+L1);

        x2_rx = y_rx_eq .* exp(-1i*omega2*n);
        x2_rx = conv(x2_rx,G2);
        x2_rx = downsample(x2_rx,M2);
        x2_est = x2_rx(S+1:S+L2);

        SINDR1(neq_idx,snr_idx) = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est-x1).^2));
        SINDR2(neq_idx,snr_idx) = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est-x2).^2));

        symerr1(neq_idx,snr_idx) = sum(qamdemod(x1_est,Q,'UnitAveragePower',true).' ~= d1);
        symerr2(neq_idx,snr_idx) = sum(qamdemod(x2_est,Q,'UnitAveragePower',true).' ~= d2);

    end
end

%% SCATTER PLOTS FOR MAIN CASES

SNR_plot = 20;
Neq_plot = 10;

snr_idx = find(SNR_test == SNR_plot);
neq_idx = find(Neq_test == Neq_plot);

SNR = SNR_test(snr_idx);
Neq = Neq_test(neq_idx);

y = upsample(y_tx,M);
y = conv(y,H_lp);

m = 0:length(y)-1;
y = real(y .* exp(1i*wc*m));

y = conv(y,c);
y = awgn(y,SNR,'measured');

m = 0:length(y)-1;
y = 2*y .* exp(-1i*wc*m);

y = conv(y,G_lp);
y = downsample(y,M);

y_rx = y(S+1:end);

[~,d_min,~,y_rx_eq] = fir_eq(y_tx,y_rx,Neq);

y_rx_eq = y_rx_eq(d_min+1:d_min+length(y_tx));

n = 0:length(y_rx_eq)-1;

x1_rx = y_rx_eq .* exp(-1i*omega1*n);
x1_rx = conv(x1_rx,G1);
x1_rx = downsample(x1_rx,M1);
x1_eq = x1_rx(S+1:S+L1);

x2_rx = y_rx_eq .* exp(-1i*omega2*n);
x2_rx = conv(x2_rx,G2);
x2_rx = downsample(x2_rx,M2);
x2_eq = x2_rx(S+1:S+L2);

figure;
subplot(3,2,1);
scatter(real(x1_ideal),imag(x1_ideal),6,'filled');
grid on; axis equal;
title('x1 ideal zero-IF');

subplot(3,2,2);
scatter(real(x2_ideal),imag(x2_ideal),6,'filled');
grid on; axis equal;
title('x2 ideal zero-IF');

subplot(3,2,3);
scatter(real(x1_channel),imag(x1_channel),6,'filled');
grid on; axis equal;
title('x1 channel without equalizer');

subplot(3,2,4);
scatter(real(x2_channel),imag(x2_channel),6,'filled');
grid on; axis equal;
title('x2 channel without equalizer');

subplot(3,2,5);
scatter(real(x1_eq),imag(x1_eq),6,'filled');
grid on; axis equal;
title('x1 equalized, SNR 20 dB, N_{eq}=10');

subplot(3,2,6);
scatter(real(x2_eq),imag(x2_eq),6,'filled');
grid on; axis equal;
title('x2 equalized, SNR 20 dB, N_{eq}=10');

%% RESULT PLOTS

figure;
subplot(2,1,1);
plot(Neq_test,SINDR1(:,snr_idx),'LineWidth',1.4); hold on;
plot(Neq_test,SINDR2(:,snr_idx),'LineWidth',1.4);
grid on;
xlabel('N_{eq}');
ylabel('SINDR [dB]');
legend('x1','x2','Location','best');
title('SINDR versus equalizer order, SNR = 20 dB');

subplot(2,1,2);
semilogy(Neq_test,symerr1(:,snr_idx)+1,'LineWidth',1.4); hold on;
semilogy(Neq_test,symerr2(:,snr_idx)+1,'LineWidth',1.4);
grid on;
xlabel('N_{eq}');
ylabel('Symbol errors + 1');
legend('x1','x2','Location','best');
title('Symbol errors versus equalizer order, SNR = 20 dB');

figure;
subplot(2,1,1);
plot(SNR_test,SINDR1(neq_idx,:),'LineWidth',1.4); hold on;
plot(SNR_test,SINDR2(neq_idx,:),'LineWidth',1.4);
grid on;
xlabel('SNR [dB]');
ylabel('SINDR [dB]');
legend('x1','x2','Location','best');
title('SINDR versus SNR, N_{eq} = 10');

subplot(2,1,2);
semilogy(SNR_test,symerr1(neq_idx,:)+1,'LineWidth',1.4); hold on;
semilogy(SNR_test,symerr2(neq_idx,:)+1,'LineWidth',1.4);
grid on;
xlabel('SNR [dB]');
ylabel('Symbol errors + 1');
legend('x1','x2','Location','best');
title('Symbol errors versus SNR, N_{eq} = 10');

%% PRINT SUMMARY

fprintf("Equalized case, SNR = 20 dB, N_eq = 10\n");
fprintf("SINDR1 = %.2f dB, symbol errors x1 = %d\n",SINDR1(neq_idx,snr_idx),symerr1(neq_idx,snr_idx));
fprintf("SINDR2 = %.2f dB, symbol errors x2 = %d\n",SINDR2(neq_idx,snr_idx),symerr2(neq_idx,snr_idx));

%% EQUALIZER FUNCTION
function [h_eq, d_min, error_min, y_rx_eq] = fir_eq(y_tx, y_rx, N_eq)

    % Length of transmitted reference signal
    L = length(y_tx);

    % Force row-vector format
    y_tx = y_tx(:).';
    y_rx = y_rx(:).';

    % Initialize outputs
    d_min = -1;
    error_min = Inf;

    % Equalizer coefficients
    h_eq = zeros(N_eq+1,1);

    % Equalized output signal
    y_rx_eq = zeros(1,L+N_eq);

    % Use same number of received samples as transmitted samples
    y_rx_use = y_rx(1:L);

    % Build Toeplitz matrix for FIR filtering
    %
    % Each row corresponds to:
    % [y_rx(n) y_rx(n-1) ... y_rx(n-N_eq)]
    %
    % Multiplying:
    % A*h_eq
    %
    % is equivalent to FIR filtering y_rx with h_eq
    col = [y_rx_use zeros(1,N_eq)];
    row = [y_rx_use(1) zeros(1,N_eq)];

    A = toeplitz(col,row);

    % Precompute least-squares matrix
    %
    % h_eq = (A'*A)^(-1)A'*y_tx
    %
    % This matrix maps the desired signal directly
    % to the least-squares equalizer coefficients
    B = (A'*A)\A';

    % Test all delays between 0 and N_eq
    for d = 0:N_eq

        % Desired delayed transmitted signal
        %
        % The equalizer output should approximate:
        %
        % y_tx(n-d)
        %
        y_tx_d = [zeros(1,d) y_tx zeros(1,N_eq-d)];

        % Least-squares equalizer coefficients
        %
        % The coefficients are complex-valued and can:
        % - rotate the constellation
        % - scale the signal
        % - cancel delayed echoes
        %
        h_eq_tmp = B*y_tx_d.';

        % Equalized received signal
        %
        % Equivalent to:
        % y_rx_eq = filter(h_eq_tmp,1,y_rx)
        %
        y_rx_eq_tmp = (A*h_eq_tmp).';

        % Mean squared error between equalized signal
        % and delayed transmitted reference signal
        err = sum(abs(y_rx_eq_tmp - y_tx_d).^2)/length(y_tx_d);

        % Store best equalizer and delay
        if err < error_min

            h_eq = h_eq_tmp;

            d_min = d;

            error_min = err;

            y_rx_eq = y_rx_eq_tmp;

        end
    end
end