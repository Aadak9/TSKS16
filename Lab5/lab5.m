clear; clc; close all;

%% PARAMETERS
% Lengths of the two transmitted 64-QAM symbol sequences
L1 = 2^13;
L2 = 2^12;

% Modulation order: 64-QAM
Q = 64;

% Sampling frequency, carrier frequency, and interpolation factor
fs = 30e6;
fc = 350e6;
M = 100;

% Upsampling factors for user 1 and user 2
M1 = 4;
M2 = 8;

% Square-root raised cosine filter parameters
rolloff = 1/3;
S = 10;

% Digital frequency shifts for the two users
omega1 = -pi/3;
omega2 = pi/6;

% Normalized carrier frequency for the zero-IF model
wc = 2*pi*fc/(M*fs);

% SNR and equalizer orders to test
SNR_test = 5:20;
Neq_test = 0:10;

% Ideal channel and multipath channel
c_ideal = 1;
c = 0.25*exp(1i*0.1*pi)*[1 zeros(1,15) 2.4 zeros(1,15) 1];

%% SIGNALS
% Generate random data symbols for both users
d1 = randi([0 Q-1],L1,1);
d2 = randi([0 Q-1],L2,1);

% Map data symbols to 64-QAM constellation points
x1 = qammod(d1,Q,'UnitAveragePower',true).';
x2 = qammod(d2,Q,'UnitAveragePower',true).';

%% FILTERS
% Pulse-shaping filters for the two transmultiplexer channels
G1 = rcosdesign(rolloff,S,M1,'sqrt')/sqrt(M1);
G2 = rcosdesign(rolloff,S,M2,'sqrt')/sqrt(M2);

% Transmit filters
H1 = M1*G1;
H2 = M2*G2;

% Lowpass filters used in the zero-IF transmitter/receiver model
G_lp = rcosdesign(rolloff,S,M,'sqrt')/sqrt(M);
H_lp = M*G_lp;

%% TRANSMUX TRANSMITTER
% Upsample, pulse-shape, and frequency-shift user 1
x1_tx = upsample(x1,M1);
x1_tx = conv(x1_tx,H1);
n = 0:length(x1_tx)-1;
x1_tx = x1_tx .* exp(1i*omega1*n);

% Upsample, pulse-shape, and frequency-shift user 2
x2_tx = upsample(x2,M2);
x2_tx = conv(x2_tx,H2);
n = 0:length(x2_tx)-1;
x2_tx = x2_tx .* exp(1i*omega2*n);

% Match signal lengths before adding both users
L = max(length(x1_tx),length(x2_tx));

x1_tx = [x1_tx zeros(1,L-length(x1_tx))];
x2_tx = [x2_tx zeros(1,L-length(x2_tx))];

% Combined transmitted transmultiplexer signal
y_tx = x1_tx + x2_tx;

%% LAB 3 REFERENCE CASE
% Ideal reference case without zero-IF model, channel, or noise
y_rx_lab3 = y_tx;

n = 0:length(y_rx_lab3)-1;

% Receiver branch for user 1
x1_rx = y_rx_lab3 .* exp(-1i*omega1*n);
x1_rx = conv(x1_rx,G1);
x1_rx = downsample(x1_rx,M1);
x1_lab3 = x1_rx(S+1:S+L1);

% Receiver branch for user 2
x2_rx = y_rx_lab3 .* exp(-1i*omega2*n);
x2_rx = conv(x2_rx,G2);
x2_rx = downsample(x2_rx,M2);
x2_lab3 = x2_rx(S+1:S+L2);

% Compute SIDR for the Lab 3 reference case
SIDR1_lab3 = 10*log10(sum(abs(x1).^2)/sum(abs(x1_lab3-x1).^2));
SIDR2_lab3 = 10*log10(sum(abs(x2).^2)/sum(abs(x2_lab3-x2).^2));

%% IDEAL ZERO-IF CASE
% Upsample the transmitted signal for RF simulation
y = upsample(y_tx,M);
y = conv(y,H_lp);

% Upconvert to RF and take the real part to model a real passband signal
m = 0:length(y)-1;
y = real(y .* exp(1i*wc*m));

% Ideal channel
y = conv(y,c_ideal);

% Downconvert from RF back to baseband
m = 0:length(y)-1;
y = 2*y .* exp(-1i*wc*m);

% Lowpass filter and downsample back to original rate
y = conv(y,G_lp);
y = downsample(y,M);

% Remove known filter delay
y_rx_ideal = y(S+1:S+length(y_tx));

n = 0:length(y_rx_ideal)-1;

% Receiver branch for user 1
x1_rx = y_rx_ideal .* exp(-1i*omega1*n);
x1_rx = conv(x1_rx,G1);
x1_rx = downsample(x1_rx,M1);
x1_ideal = x1_rx(S+1:S+L1);

% Receiver branch for user 2
x2_rx = y_rx_ideal .* exp(-1i*omega2*n);
x2_rx = conv(x2_rx,G2);
x2_rx = downsample(x2_rx,M2);
x2_ideal = x2_rx(S+1:S+L2);

% Compute SIDR and symbol errors for ideal zero-IF case
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
% Plot spectra of transmitted signal and received ideal zero-IF signal
NFFT = 2^nextpow2(length(y_tx));

YTX = fftshift(fft(y_tx,NFFT));
YRX = fftshift(fft(y_rx_ideal,NFFT));

f = linspace(-fs/2,fs/2,NFFT);

figure;

subplot(2,1,1);
%plot(f/1e6,20*log10(abs(YTX)/max(abs(YTX)) + eps));
plot(f/1e6,20*log10(abs(YTX)/max(abs(YTX))));
grid on;
xlabel('Frequency [MHz]');
ylabel('Magnitude [dB]');
title('Spectrum of y_{tx}');

subplot(2,1,2);
%plot(f/1e6,20*log10(abs(YRX)/max(abs(YRX)) + eps));
plot(f/1e6,20*log10(abs(YRX)/max(abs(YRX))));
grid on;
xlabel('Frequency [MHz]');
ylabel('Magnitude [dB]');
title('Spectrum of y_{rx} after ideal zero-IF');

%% CHANNEL WITHOUT EQUALIZER
% Repeat zero-IF transmission but now with the multipath channel
y = upsample(y_tx,M);
y = conv(y,H_lp);

m = 0:length(y)-1;
y = real(y .* exp(1i*wc*m));

% Apply multipath channel
y = conv(y,c);

m = 0:length(y)-1;
y = 2*y .* exp(-1i*wc*m);

y = conv(y,G_lp);
y = downsample(y,M);

% Remove known filter delay
y_rx_channel = y(S+1:S+length(y_tx));

n = 0:length(y_rx_channel)-1;

% Receiver branch for user 1
x1_rx = y_rx_channel .* exp(-1i*omega1*n);
x1_rx = conv(x1_rx,G1);
x1_rx = downsample(x1_rx,M1);
x1_channel = x1_rx(S+1:S+L1);

% Receiver branch for user 2
x2_rx = y_rx_channel .* exp(-1i*omega2*n);
x2_rx = conv(x2_rx,G2);
x2_rx = downsample(x2_rx,M2);
x2_channel = x2_rx(S+1:S+L2);

% Compute SINDR and symbol errors without equalization
SINDR1_channel = 10*log10(sum(abs(x1).^2)/sum(abs(x1_channel-x1).^2));
SINDR2_channel = 10*log10(sum(abs(x2).^2)/sum(abs(x2_channel-x2).^2));

err1_channel = sum(qamdemod(x1_channel,Q,'UnitAveragePower',true).' ~= d1);
err2_channel = sum(qamdemod(x2_channel,Q,'UnitAveragePower',true).' ~= d2);

fprintf("Channel without equalizer\n");
fprintf("SINDR1 = %.2f dB, symbol errors x1 = %d\n",SINDR1_channel,err1_channel);
fprintf("SINDR2 = %.2f dB, symbol errors x2 = %d\n\n",SINDR2_channel,err2_channel);

%% EQUALIZATION FOR ALL SNR AND Neq VALUES
% Allocate result matrices
SINDR1 = zeros(length(Neq_test),length(SNR_test));
SINDR2 = zeros(length(Neq_test),length(SNR_test));

symerr1 = zeros(length(Neq_test),length(SNR_test));
symerr2 = zeros(length(Neq_test),length(SNR_test));

% Loop over all tested SNR values
for snr_idx = 1:length(SNR_test)

    SNR = SNR_test(snr_idx);

    % Zero-IF transmitter
    y = upsample(y_tx,M);
    y = conv(y,H_lp);

    m = 0:length(y)-1;
    y = real(y .* exp(1i*wc*m));

    % Multipath channel and AWGN
    y = conv(y,c);
    y = awgn(y,SNR,'measured');

    % Zero-IF receiver
    m = 0:length(y)-1;
    y = 2*y .* exp(-1i*wc*m);

    y = conv(y,G_lp);
    y = downsample(y,M);

    % Remove known filter delay
    y_rx = y(S+1:end);

    % Loop over all tested equalizer orders
    for neq_idx = 1:length(Neq_test)

        Neq = Neq_test(neq_idx);

        % Design and apply FIR equalizer
        [~,d_min,~,y_rx_eq] = fir_eq(y_tx,y_rx,Neq);

        % Compensate for equalizer delay
        y_rx_eq = y_rx_eq(d_min+1:d_min+length(y_tx));

        n = 0:length(y_rx_eq)-1;

        % Receiver branch for user 1
        x1_rx = y_rx_eq .* exp(-1i*omega1*n);
        x1_rx = conv(x1_rx,G1);
        x1_rx = downsample(x1_rx,M1);
        x1_est = x1_rx(S+1:S+L1);

        % Receiver branch for user 2
        x2_rx = y_rx_eq .* exp(-1i*omega2*n);
        x2_rx = conv(x2_rx,G2);
        x2_rx = downsample(x2_rx,M2);
        x2_est = x2_rx(S+1:S+L2);

        % Compute SINDR for both users
        SINDR1(neq_idx,snr_idx) = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est-x1).^2));
        SINDR2(neq_idx,snr_idx) = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est-x2).^2));

        % Count symbol errors for both users
        symerr1(neq_idx,snr_idx) = sum(qamdemod(x1_est,Q,'UnitAveragePower',true).' ~= d1);
        symerr2(neq_idx,snr_idx) = sum(qamdemod(x2_est,Q,'UnitAveragePower',true).' ~= d2);

    end
end

%% SCATTER PLOTS FOR MAIN CASES
% Select one SNR and equalizer order for constellation plots
SNR_plot = 20;
Neq_plot = 1;

snr_idx = find(SNR_test == SNR_plot);
neq_idx = find(Neq_test == Neq_plot);

SNR = SNR_test(snr_idx);
Neq = Neq_test(neq_idx);

% Generate received signal for selected SNR
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

% Equalize received signal
[~,d_min,~,y_rx_eq] = fir_eq(y_tx,y_rx,Neq);

y_rx_eq = y_rx_eq(d_min+1:d_min+length(y_tx));

n = 0:length(y_rx_eq)-1;

% Recover user 1 after equalization
x1_rx = y_rx_eq .* exp(-1i*omega1*n);
x1_rx = conv(x1_rx,G1);
x1_rx = downsample(x1_rx,M1);
x1_eq = x1_rx(S+1:S+L1);

% Recover user 2 after equalization
x2_rx = y_rx_eq .* exp(-1i*omega2*n);
x2_rx = conv(x2_rx,G2);
x2_rx = downsample(x2_rx,M2);
x2_eq = x2_rx(S+1:S+L2);

% Plot constellations for ideal, channel-only, and equalized cases
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
% Plot SINDR versus equalizer order for SNR = 20 dB
figure;
plot(Neq_test,SINDR1(:,snr_idx),'LineWidth',1.4); hold on;
plot(Neq_test,SINDR2(:,snr_idx),'LineWidth',1.4);
grid on;
xlabel('N_{eq}');
ylabel('SINDR [dB]');
legend('x1','x2','Location','best');
title('SINDR versus equalizer order, SNR = 20 dB');

%% PRINT SUMMARY
% Print final equalized result for selected SNR and equalizer order
fprintf("Equalized case, SNR = 20 dB, N_eq = 1\n");
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

    % Initialize equalizer coefficients
    h_eq = zeros(N_eq+1,1);

    % Initialize equalized output signal
    y_rx_eq = zeros(1,L+N_eq);

    % Use same number of received samples as transmitted samples
    y_rx_use = y_rx(1:L);

    % Build Toeplitz matrix for FIR filtering
    col = [y_rx_use zeros(1,N_eq)];
    row = [y_rx_use(1) zeros(1,N_eq)];

    A = toeplitz(col,row);

    % Precompute least-squares matrix
    B = (A'*A)\A';

    % Test all possible equalizer delays
    for d = 0:N_eq

        % Desired transmitted signal with delay d
        y_tx_d = [zeros(1,d) y_tx zeros(1,N_eq-d)];

        % Compute least-squares FIR equalizer coefficients
        h_eq_tmp = B*y_tx_d.';

        % Apply equalizer 
        y_rx_eq_tmp = (A*h_eq_tmp).';

        % Compute mean squared error for this delay
        err = sum(abs(y_rx_eq_tmp - y_tx_d).^2)/length(y_tx_d);

        % Keep the equalizer delay that gives the smallest error
        if err < error_min

            h_eq = h_eq_tmp;

            d_min = d;

            error_min = err;

            y_rx_eq = y_rx_eq_tmp;

        end
    end
end