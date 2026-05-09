
%
%% TASK 6.1 CFO, PO without eachother
clc;
clear;
close all;

% Parameters
l_vals = 9:14;
L_vals = 2.^l_vals;

w0T  = 5e-5*pi;
alpha = 0.1*pi;
Q = 16;

SDR_CFO = zeros(size(L_vals));
SDR_PO  = zeros(size(L_vals));

for k = 1:length(L_vals)

    L = L_vals(k);

    % Generate random 16-QAM symbols
    data = randi([0 Q-1], L, 1);

    % 16-QAM modulation
    x = qammod(data, Q, 'UnitAveragePower', true);

    n = 0:L-1;

    %% CFO only

    y_cfo = x .* exp(1j*w0T*n.');

    SDR_CFO(k) = 10*log10(sum(abs(x).^2) / sum(abs(y_cfo - x).^2));

    figure;

    subplot(1,2,1)

    scatter(real(x), imag(x), 40, 'o');
    hold on;

    scatter(real(y_cfo), imag(y_cfo), 15, 'filled');

    axis square;
    grid on;

    title(['CFO only, L = ' num2str(L)]);
    xlabel('In-Phase');
    ylabel('Quadrature');

    legend('Original signal', 'Distorted signal');

    subplot(1,2,2)

    plot(n, real(x), 'LineWidth', 1.5);
    hold on;

    plot(n, real(y_cfo), '--', 'LineWidth', 1.2);

    grid on;

    title('Real part of signal');
    xlabel('n');
    ylabel('Amplitude');

    legend('real(x(n))', 'real(y_{cfo}(n))');

    %% Phase offset only

    y_po = x .* exp(1j*alpha);

    SDR_PO(k) = 10*log10(sum(abs(x).^2) / sum(abs(y_po - x).^2));

    figure;

    subplot(1,2,1)

    scatter(real(x), imag(x), 40, 'o');
    hold on;

    scatter(real(y_po), imag(y_po), 15, 'filled');

    axis square;
    grid on;

    title(['PO only, L = ' num2str(L)]);
    xlabel('In-Phase');
    ylabel('Quadrature');

    legend('Original signal', 'Distorted signal');

    subplot(1,2,2)

    plot(n, real(x), 'LineWidth', 1.5);
    hold on;

    plot(n, real(y_po), '--', 'LineWidth', 1.2);

    grid on;

    title('Real part of signal');
    xlabel('n');
    ylabel('Amplitude');

    legend('real(x(n))', 'real(y_{po}(n))');

end

fprintf('\nSDR RESULTS\n');

for k = 1:length(L_vals)

    fprintf('L = %6d | SDR (CFO) = %8.3f dB | SDR (PO) = %8.3f dB\n', ...
        L_vals(k), SDR_CFO(k), SDR_PO(k));

end


%SDR for PO constant since the error between x(n) and y(n) does not grow
%with time
figure;

semilogx(L_vals, SDR_CFO, '-o', 'LineWidth', 1.5);
hold on;

semilogx(L_vals, SDR_PO, '-s', 'LineWidth', 1.5);

grid on;

xlabel('Signal Length L');
ylabel('SDR (dB)');

title('SDR versus Signal Length');

legend('CFO only', 'PO only');

%}

%% TASK 6.2 Both CFO and PO present
clc;
clear;
close all;

Q = 64;

L = 2^12;

w0T = 5e-5*pi;
alpha = 0.1*pi;

SNR_dB = 30;

N_vals = [8 16 32 64 128 256];

MC = 100;

Neq = 10;

SNDR_avg = zeros(size(N_vals));

for k = 1:length(N_vals)

    N = N_vals(k);

    SNDR_temp = zeros(1, MC);

    for mc = 1:MC

        x_payload = qammod(randi([0 Q-1], L, 1), Q, 'UnitAveragePower', true);

        x = [x_payload(1:N); x_payload(1:N); x_payload];

        len = length(x);

        n = (0:len-1).';

        y = x .* exp(1j*(w0T*n + alpha));

        y = awgn(y, SNR_dB, 'measured');

        cfo_est = angle(sum(conj(y(1:N)) .* y(N+1:2*N))) / N;

        y_cfo = y .* exp(-1j*cfo_est*n);

        alpha_est = angle(sum(conj(x_payload) .* y_cfo(2*N+1:2*N+L))) / L;

        y_comp = y_cfo .* exp(-1j*alpha_est);

        [h_eq, d_min, error_min] = fir_eq(x, y_comp, Neq);

        y_eq_full = conv(y_comp, h_eq);

        y_eq = y_eq_full(d_min+1:d_min+length(x));

        idx = 2*N+1:2*N+L;

        SNDR_temp(mc) = 10*log10(sum(abs(x_payload).^2) / sum(abs(y_eq(idx) - x_payload).^2));

    end

    SNDR_avg(k) = mean(SNDR_temp);

end

figure;

semilogx(N_vals, SNDR_avg, '-o', 'LineWidth', 1.5);

grid on;

xlabel('N');

ylabel('SNDR (dB)');

title('SNDR vs N (CFO + PO + Equalizer)');

fprintf('\nRESULTS:\n');

for k = 1:length(N_vals)

    fprintf('N = %4d | SNDR = %.2f dB\n', N_vals(k), SNDR_avg(k));

end

%% TASK 6.3 CFO and PO in the receiver
clc;
clear;
close all;


%% Parameters

A1 = 1;
A2 = 1;

L1 = 2^13;
L2 = 2^12;

M = 100;

fs = 30e6;
fc = 350e6;

rolloff = 1/3;
S = 5;

Q = 64;

M1 = 4;
M2 = 8;

omega1 = -pi/3;
omega2 = pi/6;

wc = 2*pi*fc;

w0T = 5e-5*pi;
w0TH = 5e-7*pi;

alpha = 0.1*pi;

SNR_dB = 20;

Neq = 24;

%N = 128;
N = 512;

%% Channel

c = 0.25*exp(1i*0.1*pi)*[1 zeros(1,15) 2.4 zeros(1,15) 1];

%% Generate QAM symbols

x1 = qammod(randi([0 Q-1],L1,1).',Q);
x2 = qammod(randi([0 Q-1],L2,1).',Q);

%% Pulse shaping filters

G1 = rcosdesign(rolloff,S,M1,"sqrt")/sqrt(M1);
G2 = rcosdesign(rolloff,S,M2,"sqrt")/sqrt(M2);

H1 = M1*G1;
H2 = M2*G2;

%% Transmitter

x1_tx = upsample(x1,M1);
x1_tx = conv(x1_tx,H1);

n1 = 0:length(x1_tx)-1;

x1_tx = x1_tx .* exp(1i*omega1*n1);

x2_tx = upsample(x2,M2);
x2_tx = conv(x2_tx,H2);

n2 = 0:length(x2_tx)-1;

x2_tx = x2_tx .* exp(1i*omega2*n2);

%% Match lengths

L = max(length(x1_tx),length(x2_tx));

x1_tx = [x1_tx zeros(1,L-length(x1_tx))];
x2_tx = [x2_tx zeros(1,L-length(x2_tx))];

%% Combined transmitted signal

y_tx = x1_tx + x2_tx;

%% Add redundant samples for CFO estimation

y_tx = [y_tx(1:N) y_tx];

%% RF upconversion

G_lp = rcosdesign(rolloff,S,M,"sqrt")/sqrt(M);

H_lp = M*G_lp;

y_rx = upsample(y_tx,M);

y_rx = conv(y_rx,H_lp);

m = 0:length(y_rx)-1;

y_rx = y_rx .* (sqrt(2)*exp(1i*(wc/(M*fs))*m));

%% Real passband signal

y_rx = real(y_rx);

%% Channel

y_rx = conv(y_rx,c);

%% AWGN

y_rx = awgn(y_rx,SNR_dB,'measured');

%% Receiver downconversion with CFO and PO

m = 0:length(y_rx)-1;

y_rx = y_rx .* (sqrt(2)*exp(-1i*((wc/(M*fs)-w0TH)*m + alpha)));

%% Matched filter

y_rx = conv(y_rx,G_lp);

%% Downsample

y_rx = downsample(y_rx,M);

%% Remove filter transients

y_rx = y_rx(S+1:end-S);

%% CFO estimation

w0T_est = angle(sum(conj(y_rx(1:N)) .* y_rx(N+1:2*N))) / N;

fprintf('\nTrue CFO = %.6e\n',w0TH);
fprintf('Estimated CFO = %.6e\n',w0T_est);

%% CFO compensation

n = 0:length(y_rx)-1;

y_rx = y_rx .* exp(-1i*w0T_est*n);

%% Remove redundant samples

y_rx = y_rx(N+1:end);

%% Equalization

Lref = min(length(y_rx), length(y_tx)-N);

y_tx_ref = y_tx(N+1:N+Lref);

y_rx = y_rx(1:Lref);

[h_eq,d_min,error_min,y_eq] = fir_eq(y_tx_ref,y_rx,Neq);

Lcomp = min(length(y_eq)-d_min,length(y_tx_ref));

y_comp = y_eq(d_min+1:d_min+Lcomp);

%% User 1 receiver

n = 0:length(y_comp)-1;

x1_rx = y_comp .* exp(-1i*omega1*n);

x1_rx = conv(x1_rx,G1);

x1_rx = downsample(x1_rx,M1) .* (1/A1);

x1_est = x1_rx(S+1:S+L1);

%% User 2 receiver

x2_rx = y_comp .* exp(-1i*omega2*n);

x2_rx = conv(x2_rx,G2);

x2_rx = downsample(x2_rx,M2) .* (1/A2);

x2_est = x2_rx(S+1:S+L2);

%% SINDR

x1_SINDR = 10*log10(sum(abs(x1).^2) / sum(abs(x1_est-x1).^2));

x2_SINDR = 10*log10(sum(abs(x2).^2) / sum(abs(x2_est-x2).^2));

%% Symbol errors

x1_demod = qamdemod(x1,Q);
x2_demod = qamdemod(x2,Q);

x1_est_demod = qamdemod(x1_est,Q);
x2_est_demod = qamdemod(x2_est,Q);

x1_errors = sum(x1_demod ~= x1_est_demod);
x2_errors = sum(x2_demod ~= x2_est_demod);

%% Results

fprintf('\nResults for N = %d\n',N);

fprintf('Equalizer order = %d\n',Neq);

fprintf('\nUser 1\n');

fprintf('SINDR = %.2f dB\n',x1_SINDR);

fprintf('Symbol errors = %d / %d\n',x1_errors,L1);

fprintf('\nUser 2\n');

fprintf('SINDR = %.2f dB\n',x2_SINDR);

fprintf('Symbol errors = %d / %d\n',x2_errors,L2);

%% Scatter plots

Ns = 3000;

figure;

subplot(1,2,1)

hold on;
grid on;
axis square;

scatter(real(x1_est(1:Ns)), imag(x1_est(1:Ns)), 10, 'filled');
scatter(real(x1(1:Ns)), imag(x1(1:Ns)), 10, 'filled');

title('User 1: x1 vs x1\_est');
xlabel('In-Phase');
ylabel('Quadrature');
legend('x1','x1\_est');

subplot(1,2,2)

hold on;
grid on;
axis square;

scatter(real(x2_est(1:Ns)), imag(x2_est(1:Ns)), 10, 'filled');
scatter(real(x2(1:Ns)), imag(x2(1:Ns)), 10, 'filled');

title('User 2: x2 vs x2\_est');
xlabel('In-Phase');
ylabel('Quadrature');
legend('x2','x2\_est');
%}
%% Equalizer function

function [h_eq,d_min,error_min,y_rx_eq] = fir_eq(y_tx,y_rx,Neq)

L = length(y_tx);

d_min = 0;

error_min = Inf;

col = [y_rx(:); zeros(Neq,1)];

row = [y_rx(1) zeros(1,Neq)];

A = toeplitz(col,row);

B = (A'*A)\A';

for d = 0:Neq

    y_tx = y_tx(:).';
    y_rx = y_rx(:).';
    y_des = [zeros(1,d) y_tx zeros(1,Neq-d)];

    y_des = y_des(:);

    h_tmp = B * y_des;

    y_eq_tmp = A * h_tmp;

    err = mean(abs(y_eq_tmp - y_des).^2);

    if err < error_min

        error_min = err;

        h_eq = h_tmp;

        d_min = d;

        y_rx_eq = y_eq_tmp;
    end
end
y_rx_eq = y_rx_eq.';
end