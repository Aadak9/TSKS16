clear; close all; clc;

%% PARAMETERS

L1 = 2^13;
L2 = 2^12;

M  = 100;
fs = 30e6;
fc = 350e6;

rolloff = 1/3;
S  = 10;

Q  = 64;

M1 = 4;
M2 = 8;

omega1 = -pi/3;
omega2 =  pi/6;

wc = 2*pi*fc/(M*fs);

%% SIGNAL GENERATION

x1 = qammod(randi([0 Q-1], L1, 1), Q).';
x2 = qammod(randi([0 Q-1], L2, 1), Q).';

%% FILTERS

G1  = rcosdesign(rolloff, S, M1, 'sqrt') / sqrt(M1);
G2  = rcosdesign(rolloff, S, M2, 'sqrt') / sqrt(M2);

H1  = M1 * G1;
H2  = M2 * G2;

Glp = rcosdesign(rolloff, S, M, 'sqrt') / sqrt(M);
Hlp = M * Glp;

%% TRANSMITTER

x1_up = upsample(x1, M1);
x1_filt = conv(x1_up, H1);
n1 = 0:length(x1_filt)-1;
x1_mod = x1_filt .* exp(1j * omega1 * n1);

x2_up = upsample(x2, M2);
x2_filt = conv(x2_up, H2);
n2 = 0:length(x2_filt)-1;
x2_mod = x2_filt .* exp(1j * omega2 * n2);

N = max(length(x1_mod), length(x2_mod));
x1_mod = [x1_mod zeros(1, N-length(x1_mod))];
x2_mod = [x2_mod zeros(1, N-length(x2_mod))];

y_tx = x1_mod + x2_mod;

%% ZERO-IF SYSTEM

y1 = conv(y_tx, Hlp);

n = 0:length(y1)-1;
y2 = y1 .* exp(1j * wc * n);

y3 = real(y2);

%% CHANNEL

c = 0.25*exp(1j*0.1*pi)*[1 zeros(1,15) 2.4 zeros(1,15) 1];
y4 = conv(y3, c);

%% RECEIVER FRONTEND

y5 = conv(y4, Glp);

n = 0:length(y5)-1;
y_rx = y5 .* exp(-1j * wc * n);

%% DELAY COMPENSATION

delay = S*M;
y_rx = y_rx(delay+1:end);

L_eq = min(length(y_tx), length(y_rx));
y_tx_eq = y_tx(1:L_eq);
y_rx_eq = y_rx(1:L_eq);

%% EQUALIZER SWEEP

orders = 0:10;
SINDR1_vec = zeros(size(orders));
SINDR2_vec = zeros(size(orders));

for k = 1:length(orders)

    Neq = orders(k);

    [h_eq, d_min, ~] = fir_eq(y_tx_eq, y_rx_eq, Neq);

    y_rxeq = conv(y_rx, h_eq.');
    y_rxeq = y_rxeq(d_min+1 : d_min + length(y_rx));

    %% RECEIVER

    n = 0:length(y_rxeq)-1;

    y1_rx = y_rxeq .* exp(-1j * omega1 * n);
    y2_rx = y_rxeq .* exp(-1j * omega2 * n);

    z1 = conv(y1_rx, G1);
    z2 = conv(y2_rx, G2);

    z1 = z1(S*M1/2+1 : M1 : end);
    z2 = z2(S*M2/2+1 : M2 : end);

    guard = S;

    x1_est = z1(guard+1 : guard+L1);
    x2_est = z2(guard+1 : guard+L2);

    SINDR1_vec(k) = 10*log10(sum(abs(x1).^2) / sum(abs(x1_est - x1).^2));
    SINDR2_vec(k) = 10*log10(sum(abs(x2).^2) / sum(abs(x2_est - x2).^2));
end

%% RESULTS

figure;
plot(orders, SINDR1_vec, '-o', orders, SINDR2_vec, '-x');
xlabel('Equalizer order');
ylabel('SINDR (dB)');
legend('Channel 1','Channel 2');
title('SINDR vs Equalizer Order');
grid on;

%% BEST CASE CONSTELLATIONS

[h_eq, d_min, ~] = fir_eq(y_tx_eq, y_rx_eq, 10);

y_rxeq = conv(y_rx, h_eq.');
y_rxeq = y_rxeq(d_min+1 : d_min + length(y_rx));

n = 0:length(y_rxeq)-1;

y1_rx = y_rxeq .* exp(-1j * omega1 * n);
y2_rx = y_rxeq .* exp(-1j * omega2 * n);

z1 = conv(y1_rx, G1);
z2 = conv(y2_rx, G2);

z1 = z1(S*M1/2+1 : M1 : end);
z2 = z2(S*M2/2+1 : M2 : end);

x1_est = z1(guard+1 : guard+L1);
x2_est = z2(guard+1 : guard+L2);

figure;
scatterplot(x1_est);
title('x1 after equalization');

figure;
scatterplot(x2_est);
title('x2 after equalization');

%% EQUALIZER FUNCTION (UNCHANGED)

function [h_eq, d_min, error_min] = fir_eq(y_tx, y_rx, Neq)

L = length(y_tx);

d_min = 0;
error_min = Inf;

for d = 0:Neq

    y_des = [zeros(1,d), y_tx(1:L-d)];
    y_des = y_des(:);

    col = y_rx(:);
    row = [y_rx(1) zeros(1,Neq)];

    A = toeplitz(col, row);
    A = A(1:L, :);

    h_tmp = (A' * A) \ (A' * y_des);

    y_eq = A * h_tmp;

    err = sum(abs(y_eq - y_des).^2);

    if err < error_min
        error_min = err;
        h_eq = h_tmp;
        d_min = d;
    end
end

end


