%
%% TASK 6.1 CFO, PO without eachother
% This task studies CFO and PO separately for different signal lengths.
clc;
clear;
close all;

% Parameters
l_vals = 9:14;                 % Exponents used to generate signal lengths
L_vals = 2.^l_vals;            % Signal lengths L = 2^9, ..., 2^14

w0T  = 5e-5*pi;                % Normalized carrier frequency offset
alpha = 0.1*pi;                % Phase offset
Q = 16;                        % 16-QAM modulation order

SDR_CFO = zeros(size(L_vals)); % Stores SDR values when only CFO is present
SDR_PO  = zeros(size(L_vals)); % Stores SDR values when only PO is present

for k = 1:length(L_vals)

    L = L_vals(k);             % Current signal length

    % Generate random 16-QAM symbols
    data = randi([0 Q-1], L, 1);

    % 16-QAM modulation
    x = qammod(data, Q, 'UnitAveragePower', true);

    n = 0:L-1;                 % Sample index vector

    %% CFO only
    % Apply only carrier frequency offset.
    % The phase rotation increases with time/sample index n.

    y_cfo = x .* exp(1j*w0T*n.');

    % Compute SDR between original signal and CFO-distorted signal
    SDR_CFO(k) = 10*log10(sum(abs(x).^2) / sum(abs(y_cfo - x).^2));

    figure;

    subplot(1,2,1)

    % Plot original and CFO-distorted constellations
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

    % Plot real part before and after CFO
    plot(n, real(x), 'LineWidth', 1.5);
    hold on;

    plot(n, real(y_cfo), '--', 'LineWidth', 1.2);

    grid on;

    title('Real part of signal');
    xlabel('n');
    ylabel('Amplitude');

    legend('real(x(n))', 'real(y_{cfo}(n))');

    %% Phase offset only
    % Apply only phase offset.
    % This rotates all constellation points by the same fixed angle.

    y_po = x .* exp(1j*alpha);

    % Compute SDR between original signal and PO-distorted signal
    SDR_PO(k) = 10*log10(sum(abs(x).^2) / sum(abs(y_po - x).^2));

    figure;

    subplot(1,2,1)

    % Plot original and phase-offset constellations
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

    % Plot real part before and after phase offset
    plot(n, real(x), 'LineWidth', 1.5);
    hold on;

    plot(n, real(y_po), '--', 'LineWidth', 1.2);

    grid on;

    title('Real part of signal');
    xlabel('n');
    ylabel('Amplitude');

    legend('real(x(n))', 'real(y_{po}(n))');

end

% Print SDR values for all tested signal lengths
fprintf('\nSDR RESULTS\n');

for k = 1:length(L_vals)

    fprintf('L = %6d | SDR (CFO) = %8.3f dB | SDR (PO) = %8.3f dB\n', ...
        L_vals(k), SDR_CFO(k), SDR_PO(k));

end


% SDR for PO constant since the error between x(n) and y(n) does not grow
% with time.
% CFO SDR decreases with L because CFO rotation accumulates over time.
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
% This task estimates and compensates CFO and PO using a repeated training block.
%clc;
%clear;
%close all;

Q = 64; % 64-QAM modulation

L = 2^12; % Payload length

w0T = 5e-5*pi; % Carrier frequency offset (CFO)
alpha = 0.1*pi; % Constant phase offset (PO)

SNR_dB = 30; % AWGN signal-to-noise ratio

N_vals = [8 16 32 64 128 256]; % Repeated-block lengths for CFO estimation

MC = 100; % Number of Monte Carlo simulations

Neq = 10; % Equalizer order

SNDR_avg = zeros(size(N_vals)); % Storage for average SNDR values

for k = 1:length(N_vals)

    N = N_vals(k); % Current repeated-block length

    SNDR_temp = zeros(1, MC); % Storage for temporary SNDR values

    for mc = 1:MC

        % Generate random 64-QAM payload symbols
        x_payload = qammod(randi([0 Q-1], L, 1), Q, 'UnitAveragePower', true);

        % Create repeated training sequence:
        % first N samples repeated twice for CFO estimation
        x = [x_payload(1:N); x_payload(1:N); x_payload];

        len = length(x); % Total signal length

        n = (0:len-1).'; % Sample indices

        % Apply CFO and PO
        % CFO creates a time-varying phase rotation.
        % PO creates a constant phase rotation.
        y = x .* exp(1j*(w0T*n + alpha));

        % Add AWGN
        y = awgn(y, SNR_dB, 'measured');

        % CFO estimation using repeated blocks
        % Since the first N samples are repeated, the phase difference
        % between the two blocks can be used to estimate CFO.
        cfo_est = angle(sum(conj(y(1:N)) .* y(N+1:2*N))) / N;

        % CFO compensation
        % Multiply by the opposite phase rotation.
        y_cfo = y .* exp(-1j*cfo_est*n);

        % Phase offset estimation after CFO compensation
        % The training blocks are skipped so only the payload is used here.
        alpha_est = angle(sum(conj(x_payload) .* y_cfo(2*N+1:2*N+L))) / L;%%Make sure to skip training blocks

        % Phase offset compensation
        y_comp = y_cfo .* exp(-1j*alpha_est);

        
        % Design FIR equalizer using transmitted signal as reference
        [h_eq, d_min, error_min] = fir_eq(x, y_comp, Neq);

        % Equalize received signal
        y_eq_full = conv(y_comp, h_eq);

        % Remove equalizer delay
        y_eq = y_eq_full(d_min+1:d_min+length(x));

        % Extract payload region
        idx = 2*N+1:2*N+L;

        % Compute SNDR
        % Compares compensated received payload with original payload.
        SNDR_temp(mc) = 10*log10(sum(abs(x_payload).^2) / sum(abs(y_eq(idx) - x_payload).^2));

    end

    % Average SNDR over all Monte Carlo runs
    SNDR_avg(k) = mean(SNDR_temp);

end

figure;

% Plot average SNDR versus repeated-block length N
semilogx(N_vals, SNDR_avg, '-o', 'LineWidth', 1.5);

grid on;

xlabel('N');

ylabel('SNDR (dB)');

title('SNDR vs N (CFO + PO + Equalizer)');

fprintf('\nRESULTS:\n');

for k = 1:length(N_vals)

    % Print SNDR results for each N
    fprintf('N = %4d | SNDR = %.2f dB\n', N_vals(k), SNDR_avg(k));

end

%% TASK 6.3 CFO and PO in the receiver
% This task introduces CFO and PO in the full zero-IF receiver chain.
%clc;
%clear;
%close all;


%% Parameters

A1 = 1;                         % Gain factor for user 1
A2 = 1;                         % Gain factor for user 2

L1 = 2^13;                      % Number of symbols for user 1
L2 = 2^12;                      % Number of symbols for user 2

M = 100;                        % Interpolation/downsampling factor for zero-IF model

fs = 30e6;                      % Sampling frequency
fc = 350e6;                     % Carrier frequency

rolloff = 1/3;                  % Raised cosine rolloff factor
S = 5;                          % Filter span

Q = 64;                         % 64-QAM modulation

M1 = 4;                         % Upsampling factor for user 1
M2 = 8;                         % Upsampling factor for user 2

omega1 = -pi/3;                 % Digital frequency shift for user 1
omega2 = pi/6;                  % Digital frequency shift for user 2

wc = 2*pi*fc;                   % Angular carrier frequency

w0T = 5e-5*pi;                  % CFO at lower symbol-rate model
w0TH = 5e-7*pi;                 % CFO at higher-rate zero-IF model

alpha = 0.1*pi;                 % Phase offset

SNR_dB = 20;                    % AWGN SNR

Neq = 24;                       % Equalizer order

%N = 128;
N = 512;                        % Number of redundant samples used for CFO estimation

%% Channel

% Multipath channel with three taps at delays 0, 16, and 32
c = 0.25*exp(1i*0.1*pi)*[1 zeros(1,15) 2.4 zeros(1,15) 1];

%% Generate QAM symbols

% Generate 64-QAM symbols for both users
x1 = qammod(randi([0 Q-1],L1,1).',Q,'UnitAveragePower',true);
x2 = qammod(randi([0 Q-1],L2,1).',Q,'UnitAveragePower',true);

%% Pulse shaping filters

% Square-root raised cosine filters for pulse shaping
G1 = rcosdesign(rolloff,S,M1,"sqrt")/sqrt(M1);
G2 = rcosdesign(rolloff,S,M2,"sqrt")/sqrt(M2);

% Transmit filters
H1 = M1*G1;
H2 = M2*G2;

%% Transmitter

% User 1 transmitter branch:
% upsample, pulse-shape, and shift to frequency omega1
x1_tx = upsample(x1,M1);
x1_tx = conv(x1_tx,H1);

n1 = 0:length(x1_tx)-1;

x1_tx = x1_tx .* exp(1i*omega1*n1);

% User 2 transmitter branch:
% upsample, pulse-shape, and shift to frequency omega2
x2_tx = upsample(x2,M2);
x2_tx = conv(x2_tx,H2);

n2 = 0:length(x2_tx)-1;

x2_tx = x2_tx .* exp(1i*omega2*n2);

%% Match lengths

% Make both transmitted user signals equally long before adding them
L = max(length(x1_tx),length(x2_tx));

x1_tx = [x1_tx zeros(1,L-length(x1_tx))];
x2_tx = [x2_tx zeros(1,L-length(x2_tx))];

%% Combined transmitted signal

% Add both users to form the transmultiplexer output
y_tx = x1_tx + x2_tx;

%% Add redundant samples for CFO estimation

% Copy the first N samples and prepend them.
% These repeated samples are later used to estimate CFO.
y_tx = [y_tx(1:N) y_tx];

%% RF upconversion

% Lowpass interpolation filter for zero-IF model
G_lp = rcosdesign(rolloff,S,M,"sqrt")/sqrt(M);

H_lp = M*G_lp;

% Upsample transmitted baseband signal
y_rx = upsample(y_tx,M);

% Interpolation filtering
y_rx = conv(y_rx,H_lp);

m = 0:length(y_rx)-1;

% Upconvert complex baseband signal to RF
y_rx = y_rx .* (sqrt(2)*exp(1i*(wc/(M*fs))*m));

%% Real passband signal

% Physical RF signals are real-valued
y_rx = real(y_rx);

%% Channel

% Apply multipath channel
y_rx = conv(y_rx,c);

%% AWGN

% Add white Gaussian noise
y_rx = awgn(y_rx,SNR_dB,'measured');

%% Receiver downconversion with CFO and PO

m = 0:length(y_rx)-1;

% Downconvert from RF to baseband.
% The receiver oscillator includes CFO and phase offset.
y_rx = y_rx .* (sqrt(2)*exp(-1i*((wc/(M*fs)-w0TH)*m + alpha)));

%% Matched filter

% Matched filtering after downconversion
y_rx = conv(y_rx,G_lp);

%% Downsample

% Return to lower sampling rate
y_rx = downsample(y_rx,M);

%% Remove filter transients

% Remove transient samples caused by filtering
y_rx = y_rx(S+1:end-S);

%% CFO estimation

% Estimate CFO from the repeated N-sample block
w0T_est = angle(sum(conj(y_rx(1:N)) .* y_rx(N+1:2*N))) / N;

fprintf('\nTrue CFO = %.6e\n',w0TH);
fprintf('Estimated CFO = %.6e\n',w0T_est);

%% CFO compensation

n = 0:length(y_rx)-1;

% Remove estimated CFO from received signal
y_rx = y_rx .* exp(-1i*w0T_est*n);

%% Remove redundant samples

% Remove the N copied samples after CFO estimation
y_rx = y_rx(N+1:end);

%% Equalization

% Choose common length for transmitted reference and received signal
Lref = min(length(y_rx), length(y_tx)-N);

% Reference signal without redundant CFO-training samples
y_tx_ref = y_tx(N+1:N+Lref);

% Truncate received signal to same reference length
y_rx = y_rx(1:Lref);

% Design and apply FIR equalizer
[h_eq,d_min,error_min,y_eq] = fir_eq(y_tx_ref,y_rx,Neq);

% Determine valid comparison length after delay compensation
Lcomp = min(length(y_eq)-d_min,length(y_tx_ref));

% Remove equalizer delay
y_comp = y_eq(d_min+1:d_min+Lcomp);

%% User 1 receiver

n = 0:length(y_comp)-1;

% Shift user 1 back to baseband
x1_rx = y_comp .* exp(-1i*omega1*n);

% Matched filter for user 1
x1_rx = conv(x1_rx,G1);

% Downsample and compensate gain
x1_rx = downsample(x1_rx,M1) .* (1/A1);

% Remove filter delay and keep only transmitted symbol interval
x1_est = x1_rx(S+1:S+L1);

%% User 2 receiver

% Shift user 2 back to baseband
x2_rx = y_comp .* exp(-1i*omega2*n);

% Matched filter for user 2
x2_rx = conv(x2_rx,G2);

% Downsample and compensate gain
x2_rx = downsample(x2_rx,M2) .* (1/A2);

% Remove filter delay and keep only transmitted symbol interval
x2_est = x2_rx(S+1:S+L2);

%% SINDR

% Compute SINDR for user 1
x1_SINDR = 10*log10(sum(abs(x1).^2) / sum(abs(x1_est-x1).^2));

% Compute SINDR for user 2
x2_SINDR = 10*log10(sum(abs(x2).^2) / sum(abs(x2_est-x2).^2));

%% Symbol errors

% Demodulate original transmitted symbols
x1_demod = qamdemod(x1,Q);
x2_demod = qamdemod(x2,Q);

% Demodulate estimated received symbols
x1_est_demod = qamdemod(x1_est,Q);
x2_est_demod = qamdemod(x2_est,Q);

% Count symbol errors
x1_errors = sum(x1_demod ~= x1_est_demod);
x2_errors = sum(x2_demod ~= x2_est_demod);

%% Results

% Print final performance results
fprintf('\nResults for N = %d\n',N);

fprintf('Equalizer order = %d\n',Neq);

fprintf('\nUser 1\n');

fprintf('SINDR = %.2f dB\n',x1_SINDR);

fprintf('Symbol errors = %d / %d\n',x1_errors,L1);

fprintf('\nUser 2\n');

fprintf('SINDR = %.2f dB\n',x2_SINDR);

fprintf('Symbol errors = %d / %d\n',x2_errors,L2);

%% Scatter plots

Ns = 3000;                      % Number of symbols shown in scatter plots

figure;

subplot(1,2,1)

hold on;
grid on;
axis square;

% Compare estimated and original constellation for user 1
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

% Compare estimated and original constellation for user 2
scatter(real(x2_est(1:Ns)), imag(x2_est(1:Ns)), 10, 'filled');
scatter(real(x2(1:Ns)), imag(x2(1:Ns)), 10, 'filled');

title('User 2: x2 vs x2\_est');
xlabel('In-Phase');
ylabel('Quadrature');
legend('x2','x2\_est');
%}

%% Equalizer function
% This function designs a least-squares FIR equalizer.
% It tries different delays and selects the delay with the smallest error.

function [h_eq,d_min,error_min,y_rx_eq] = fir_eq(y_tx,y_rx,Neq)

% Length of transmitted reference signal
L = length(y_tx);

% Initialize best delay
d_min = 0;

% Initialize minimum error as infinity
error_min = Inf;

% First column of Toeplitz convolution matrix
col = [y_rx(:); zeros(Neq,1)];

% First row of Toeplitz convolution matrix
row = [y_rx(1) zeros(1,Neq)];

% Toeplitz matrix represents FIR filtering as matrix multiplication
A = toeplitz(col,row);

% Precompute least-squares pseudoinverse matrix
B = (A'*A)\A';

% Try all possible equalizer delays from 0 to Neq
for d = 0:Neq

    % Force signals to row-vector format
    y_tx = y_tx(:).';
    y_rx = y_rx(:).';

    % Desired equalizer output:
    % delayed transmitted reference signal
    y_des = [zeros(1,d) y_tx zeros(1,Neq-d)];

    % Convert desired signal to column vector
    y_des = y_des(:);

    % Compute least-squares equalizer coefficients for this delay
    h_tmp = B * y_des;

    % Apply equalizer using matrix multiplication
    y_eq_tmp = A * h_tmp;

    % Compute mean squared error for this delay
    err = mean(abs(y_eq_tmp - y_des).^2);

    % Store equalizer if this delay gives lower error
    if err < error_min

        error_min = err;

        h_eq = h_tmp;

        d_min = d;

        y_rx_eq = y_eq_tmp;
    end
end

% Return equalized signal as row vector
y_rx_eq = y_rx_eq.';
end