clear; clc; close all;

%% PARAMETERS

L1 = 2^13;
L2 = 2^12;

M = 100;
fs = 30e6;

rolloff = 1/3;
S = 5;

Q = 64;

M1 = 4;
M2 = 8;

omega1 = -pi/3;
omega2 = pi/6;

SNR_test = 5:20;

c = 0.25*exp(1i*0.1*pi)*[1 zeros(1,15) 2.4 zeros(1,15) 1];

SINDR_x1 = zeros(1,length(SNR_test));
SINDR_x2 = zeros(1,length(SNR_test));

%% SIGNALS

x1 = qammod(randi([0 Q-1], L1, 1), Q).';
x2 = qammod(randi([0 Q-1], L2, 1), Q).';

%% FILTERS

G1 = rcosdesign(rolloff, S, M1, "sqrt")/sqrt(M1);
G2 = rcosdesign(rolloff, S, M2, "sqrt")/sqrt(M2);

H1 = M1*G1;
H2 = M2*G2;

G_lp = rcosdesign(rolloff, S, M, "sqrt")/sqrt(M);
H_lp = M*G_lp;

%% TRANSMITTER

x1_tx = upsample(x1, M1);
x1_tx = conv(x1_tx, H1);
x1_tx = x1_tx .* exp(1i*omega1*(0:length(x1_tx)-1));

x2_tx = upsample(x2, M2);
x2_tx = conv(x2_tx, H2);
x2_tx = x2_tx .* exp(1i*omega2*(0:length(x2_tx)-1));

L = max(length(x1_tx), length(x2_tx));

x1_tx = [x1_tx zeros(1,L-length(x1_tx))];
x2_tx = [x2_tx zeros(1,L-length(x2_tx))];

y_tx = x1_tx + x2_tx;

y_tx_if = y_tx;   % keep for spectrum

%% TX SPECTRUM (PLOT 1)

figure;
pwelch(y_tx_if,[],[],[],fs*M);
title('Transmit Spectrum y_{tx}');
grid on;

%% ZERO-IF + CHANNEL

y = upsample(y_tx, M);
y = conv(y, H_lp);

y = y .* (sqrt(2)*exp(1i*(0:length(y)-1)*0.1));

y = real(y);

y = conv(y, c);

y_rx0 = y;

%% CHANNEL SPECTRUM (PLOT 2)

figure;
pwelch(y_rx0,[],[],[],fs*M);
title('Channel Output Spectrum');
grid on;

%% EQUALIZER + SNR LOOP

for idx = 1:length(SNR_test)

    SNR = SNR_test(idx);

    y_rx = awgn(y_rx0, SNR, "measured");

    y_rx = y_rx .* (sqrt(2)*exp(-1i*(0:length(y_rx)-1)*0.1));

    y_rx = conv(y_rx, G_lp);
    y_rx = downsample(y_rx, M);
    y_rx = y_rx(S+1:end-S);

    y_tx_ext = [y_tx zeros(1,length(y_rx)-length(y_tx))];

    [h_eq, d_min, ~, y_rx_eq] = eq(y_tx_ext, y_rx, M);

    y_rx_eq = y_rx_eq(d_min+1:d_min+length(y_tx));

    %% RECEIVER

    x1_rx = y_rx_eq .* exp(-1i*omega1*(0:length(y_rx_eq)-1));
    x1_rx = conv(x1_rx, G1);
    x1_rx = downsample(x1_rx, M1);

    x2_rx = y_rx_eq .* exp(-1i*omega2*(0:length(y_rx_eq)-1));
    x2_rx = conv(x2_rx, G2);
    x2_rx = downsample(x2_rx, M2);

    guard = S;

    x1_est = x1_rx(guard+1:guard+L1);
    x2_est = x2_rx(guard+1:guard+L2);

    SINDR_x1(idx) = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est-x1).^2));
    SINDR_x2(idx) = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est-x2).^2));

end

%% RESULTS

figure;

subplot(2,1,1);
hold on; 
grid on;
scatter(real(x1_est), imag(x1_est), 10, 'r', 'filled');
scatter(real(x1), imag(x1), 20, 'k', 'filled');
axis equal;
legend('x1','x1 est');
title('x1 Constellation');

subplot(2,1,2);
hold on; 
grid on;
scatter(real(x2_est), imag(x2_est), 10, 'r', 'filled');
scatter(real(x2), imag(x2), 20, 'k', 'filled');
axis equal;
legend('x2','x2 est');
title('x2 Constellation');

%% SINDR PLOT (PLOT 3)

figure;
plot(SNR_test, SINDR_x1); 
hold on;
plot(SNR_test, SINDR_x2);
grid on;
legend('x1','x2');
title('SINDR vs SNR');
xlabel('SNR (dB)');
ylabel('SINDR (dB)');

%% EQUALIZER FUNCTION

function [h_eq, d_min, error_min, y_rx_eq] = eq(y_tx, y_rx, N_eq)

L = length(y_tx);

d_min = 0;
error_min = Inf;

col = [y_rx zeros(1,N_eq)];
row = [y_rx(1) zeros(1,N_eq)];

A = toeplitz(col,row);
B = inv(A'*A)*A';

for d = 0:N_eq

    y_tx_d = [zeros(1,d) y_tx zeros(1,N_eq-d)];
    h_eq_tmp = B*y_tx_d.';

    y_rx_eq_tmp = (A*h_eq_tmp).';

    err = sum(abs(y_rx_eq_tmp - y_tx_d).^2);

    if err < error_min
        error_min = err;
        h_eq = h_eq_tmp;
        d_min = d;
        y_rx_eq = y_rx_eq_tmp;
    end
end

end