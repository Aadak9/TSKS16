close all;
clear;
clc;

%% Parameters

M = 100;                        % Interpolation/downsampling factor in RF model
M1 = 4;                         % Upsampling factor for signal 1
M2 = 8;                         % Upsampling factor for signal 2

L1 = 2^13;                      % Number of symbols for signal 1
L2 = 2^12;                      % Number of symbols for signal 2

fs = 30e6;                      % Baseband sampling frequency
fc = 350e6;                     % Carrier frequency
fs_RF = M*fs;                   % Sampling frequency after interpolation

rolloff = 1/3;                  % RRC rolloff factor

A1 = 0.05;                      % Amplitude scaling for signal 1
A2 = 1;                         % Amplitude scaling for signal 2

S1 = 30;                        % Filter span for signal 1
S2 = 30;                        % Filter span for signal 2
S_LP = 20;                      % Filter span for RF digital model

Q = 16;                         % 16-QAM
noise_SNR = 50;                 % AWGN SNR

c = zeros(33,1);                % Multipath channel
c(1) = 0.25;
c(17) = 2.4*0.25;
c(33) = 0.25;

wT = 5e-7*pi;                   % CFO value from Lab 7
alpha = 0.1*pi;                 % Phase offset
N_offset = 512;                 % Repeated block length for CFO estimation

g1 = 0.99;                      % I-branch gain mismatch
g2 = 1.02;                      % Q-branch gain mismatch

phi1 = 0.5;                     % I-branch phase mismatch
phi2 = 0.55;                    % Q-branch phase mismatch

omega1 = -pi/3;                 % Frequency shift for signal 1
omega2 = pi/6;                  % Frequency shift for signal 2

omega_c = 2*pi*fc/(M*fs);       % Normalized RF carrier frequency

N_eq = 6;                       % Equalizer order

rng(1);                         % Fixed random seed for repeatable results

%% Generate QAM symbols

x1 = qammod(randi([0 Q-1],L1,1),Q,'UnitAveragePower',true);
x2 = qammod(randi([0 Q-1],L2,1),Q,'UnitAveragePower',true);

%% Pulse shaping filters

G1 = rcosdesign(rolloff,S1,M1,'sqrt')/sqrt(M1);
G2 = rcosdesign(rolloff,S2,M2,'sqrt')/sqrt(M2);

H1 = M1*G1;
H2 = M2*G2;

G_LP = rcosdesign(rolloff,S_LP,M,'sqrt')/sqrt(M);
H_LP = M*G_LP;

%% Bandpass filter design

% The bandpass filter models the analog RF bandpass filter before demodulation.
% It suppresses out-of-band nonlinear distortion products.

Fstop1 = fc - 30e6;
Fpass1 = fc - 20e6;
Fpass2 = fc + 20e6;
Fstop2 = fc + 30e6;

dev = [0.001 0.001 0.001];

[n_bpf,fo_bpf,ao_bpf,w_bpf] = firpmord([Fstop1 Fpass1 Fpass2 Fstop2],[0 1 0],dev,fs_RF);

n_bpf = n_bpf + rem(n_bpf,2);

b_bpf = firpm(n_bpf,fo_bpf,ao_bpf,w_bpf);

%% Case settings

case_names = {'Case 1: Nonlinearity alone, no compensation','Case 2: Nonlinearity alone, with linearization and equalization','Case 3: Nonlinearity, IQM, CFO, PO, channel, bandpass filter, no compensation','Case 4: Nonlinearity, IQM, CFO, PO, channel, bandpass filter, with compensation'};

case_NL_ON = [1 1 1 1];

case_LIN_COMP = [0 1 0 1];

case_BPF_ON = [0 0 1 1];

case_IQM_ON = [0 0 1 1];
case_IQM_COMP = [0 0 0 1];

case_CFO_ON = [0 0 1 1];
case_CFO_COMP = [0 0 0 1];

case_CHANNEL_ON = [0 0 1 1];

case_EQ_ON = [0 1 0 1];

res = struct();

%% Run main cases

for case_index = 1:length(case_names)

    fprintf('\n%s\n',case_names{case_index});

    NL_ON = case_NL_ON(case_index);

    LIN_COMP = case_LIN_COMP(case_index);

    BPF_ON = case_BPF_ON(case_index);

    IQM_ON = case_IQM_ON(case_index);
    IQM_COMP = case_IQM_COMP(case_index);

    CFO_ON = case_CFO_ON(case_index);
    CFO_COMP = case_CFO_COMP(case_index);

    CHANNEL_ON = case_CHANNEL_ON(case_index);

    EQ_ON = case_EQ_ON(case_index);

    %% Transmitter

    % Upsample both QAM symbol sequences
    x1_up = upsample(x1,M1);
    x2_up = upsample(x2,M2);

    % Apply pulse shaping and amplitude scaling
    y1 = A1*conv(x1_up,H1);
    y2 = A2*conv(x2_up,H2);

    % Remove pulse-shaping filter delay
    y1 = y1((length(H1)-1)/2+1:end-(length(H1)-1)/2);
    y2 = y2((length(H2)-1)/2+1:end-(length(H2)-1)/2);

    % Make the two transmitted signals equally long
    if length(y1) > length(y2)
        y2(end+1:length(y1)) = 0;
    else
        y1(end+1:length(y2)) = 0;
    end

    % Frequency-shift both signals to their subbands
    N_signal = length(y1);
    n = (0:N_signal-1).';

    factor1 = exp(1i*omega1*n);
    factor2 = exp(1i*omega2*n);

    y1 = y1.*factor1;
    y2 = y2.*factor2;

    % Add both signals
    ytx = y1 + y2;

    % Add repeated block for CFO estimation when CFO is present
    if CFO_ON == 1
        ytx = [ytx(1:N_offset); ytx];
    end

    %% RF upconversion

    % Upsample to RF model sampling frequency
    ytx_up = upsample(ytx,M);

    % Interpolation filter
    ytx_filtered = conv(ytx_up,H_LP);
    ytx_filtered = ytx_filtered((length(H_LP)-1)/2+1:end-(length(H_LP)-1)/2);

    % RF-rate time index
    N_upsampled = length(ytx_filtered);
    m = (0:N_upsampled-1).';

    % Upconvert to RF
    carrier_up = sqrt(2)*exp(1i*omega_c*m);

    ytx_modulated = ytx_filtered.*carrier_up;

    % Real RF signal before nonlinearity
    r = real(ytx_modulated);

    res(case_index).r_before_nonlinearity = r;

    %% Nonlinearity

    % Polynomial nonlinearity:
    % y(m) = a0 + a1*r(m) + a2*r^2(m) + a3*r^3(m)
    if NL_ON == 1

        a0 = 0.01;
        a1 = 1;
        a2 = -0.2/max(abs(r));
        a3 = 0.15/(max(abs(r))^2);

        yrx = a0 + a1*r + a2*r.^2 + a3*r.^3;

    else

        a0 = 0;
        a1 = 1;
        a2 = 0;
        a3 = 0;

        yrx = r;

    end

    res(case_index).yrx_after_nonlinearity = yrx;

    %% Channel

    % Apply multipath channel when enabled
    if CHANNEL_ON == 1
        yrx = conv(yrx,c);
        yrx = yrx((length(c)-1)/2+1:end-(length(c)-1)/2);
    end

    %% AWGN

    % Add AWGN
    if noise_SNR ~= 0
        yrx = awgn(yrx,noise_SNR,'measured');
    end

    res(case_index).yrx_before_bpf = yrx;

    %% Bandpass filter

    % Apply RF bandpass filter before demodulation
    if BPF_ON == 1
        yrx = conv(yrx,b_bpf);
        yrx = yrx(n_bpf/2+1:end-n_bpf/2);
    end

    res(case_index).yrx_after_bpf = yrx;

    %% Receiver downconversion

    % Receiver oscillator frequency
    if CFO_ON == 1
        omega_down = omega_c - wT/M;
    else
        omega_down = omega_c;
    end

    % Receiver demodulator with or without IQ mismatch
    if IQM_ON == 1
        if CFO_ON == 1
            carrier_down = sqrt(2)*(g1*cos(omega_down*m + phi1 - alpha) - 1i*g2*sin(omega_down*m + phi2 - alpha));
        else
            carrier_down = sqrt(2)*(g1*cos(omega_down*m + phi1) - 1i*g2*sin(omega_down*m + phi2));
        end
    else
        if CFO_ON == 1
            carrier_down = sqrt(2)*exp(-1i*(omega_down*m - alpha));
        else
            carrier_down = sqrt(2)*exp(-1i*omega_down*m);
        end
    end

    % Match length after possible bandpass filtering
    Lmix = min(length(yrx),length(carrier_down));

    yrx = yrx(1:Lmix);
    carrier_down = carrier_down(1:Lmix);

    % Downconvert to complex baseband
    yrx_modulated = yrx.*carrier_down;

    %% Receiver lowpass filtering and downsampling

    % Lowpass filter after demodulation
    yrx_filtered = conv(yrx_modulated,G_LP);
    yrx_filtered = yrx_filtered((length(G_LP)-1)/2+1:end-(length(G_LP)-1)/2);

    % Downsample back to baseband rate
    yrx = downsample(yrx_filtered,M);

    res(case_index).yrx_before_linearization = yrx;

    %% Linearization

    % Baseband compensation of cubic nonlinearity:
    % yrx(n) = yrx(n) - 3*a3*|yrx(n)|^2*yrx(n)
    if LIN_COMP == 1
    yrx = yrx - a3*abs(yrx).^2.*yrx;
    end

    res(case_index).yrx_after_linearization = yrx;

    %% IQM compensation

    % First-order IQM compensation from Lab 7
    if IQM_ON == 1 && IQM_COMP == 1

        yr = real(yrx);
        yi = imag(yrx);

        Pr = mean(yr.^2);
        Pi = mean(yi.^2);

        g_est = sqrt(Pi/Pr);

        yp = -yr.*yi;
        Pyp = mean(yp);

        phi_est = Pyp/Pr;

        yr = yr*g_est;
        yi = yi + yr*phi_est;

        yrx = yr + 1i*yi;

    end

    res(case_index).yrx_after_iqm_comp = yrx;

    %% CFO estimation

    % Estimate CFO from the repeated block
    if CFO_ON == 1
        CFO_est = angle(sum(conj(yrx(1:N_offset)).*yrx(N_offset+1:2*N_offset)))/N_offset;
    else
        CFO_est = 0;
    end

    res(case_index).CFO_est = CFO_est;

    %% CFO compensation

    % Remove estimated CFO
    if CFO_ON == 1 && CFO_COMP == 1
        n_temp = (0:length(yrx)-1).';
        yrx = yrx.*exp(-1i*CFO_est*n_temp);
    end

    %% Remove repeated CFO block

    % Remove the copied block after CFO estimation
    if CFO_ON == 1
        yrx = yrx(N_offset+1:end);
        ytx = ytx(N_offset+1:end);
    end

    %% Equalization

    % Equalizer compensates channel, PO, and remaining linear distortion
    if EQ_ON == 1

        Lref = min(length(ytx),length(yrx));

        ytx_ref = ytx(1:Lref);
        yrx_ref = yrx(1:Lref);

        [h_eq,d_min,error_min,yrx_eq_full] = fir_eq(ytx_ref,yrx_ref,N_eq);

        Lcomp = min(length(yrx_eq_full)-d_min,length(ytx_ref));

        yrx_eq = yrx_eq_full(d_min+1:d_min+Lcomp).';

    else

        yrx_eq = yrx;

    end

    res(case_index).yrx_after_equalizer = yrx_eq;

    %% Receiver demultiplexing

    % Time index for receiver branch separation
    N_rx = length(yrx_eq);
    n_rx = (0:N_rx-1).';

    factor1_rx = exp(1i*omega1*n_rx);
    factor2_rx = exp(1i*omega2*n_rx);

    % Shift each signal back to baseband
    y_rx1 = yrx_eq.*conj(factor1_rx);
    y_rx2 = yrx_eq.*conj(factor2_rx);

    % Matched filtering
    r1 = conv(y_rx1,G1);
    r2 = conv(y_rx2,G2);

    % Remove matched-filter delays
    r1 = r1((length(G1)-1)/2+1:end-(length(G1)-1)/2);
    r2 = r2((length(G2)-1)/2+1:end-(length(G2)-1)/2);

    % Downsample to symbol rate
    x1_est = downsample(r1,M1);
    x2_est = downsample(r2,M2);

    % Undo amplitude scaling
    if A1 ~= 0
        x1_est = x1_est/A1;
    end

    x2_est = x2_est/A2;

    %% SINDR and symbol errors

    % Use common lengths for transmitted and estimated symbol sequences
    Lx1 = min(length(x1),length(x1_est));
    Lx2 = min(length(x2),length(x2_est));

    x1_use = x1(1:Lx1);
    x2_use = x2(1:Lx2);

    x1_est = x1_est(1:Lx1);
    x2_est = x2_est(1:Lx2);

    % Signal 1 SINDR and symbol errors
    if A1 ~= 0

        b1 = (x1_use'*x1_est)/(x1_est'*x1_est);

        x1_est_aligned = b1*x1_est;

        SINDR1 = 10*log10(sum(abs(x1_use).^2)/sum(abs(x1_use - x1_est_aligned).^2));

        x1_demod = qamdemod(x1_use,Q,'UnitAveragePower',true);

        x1_est_demod = qamdemod(x1_est_aligned,Q,'UnitAveragePower',true);

        errors1 = sum(x1_demod ~= x1_est_demod);

    else

        x1_est_aligned = x1_est;

        SINDR1 = NaN;

        errors1 = NaN;

    end

    % Signal 2 SINDR and symbol errors
    b2 = (x2_use'*x2_est)/(x2_est'*x2_est);

    x2_est_aligned = b2*x2_est;

    SINDR2 = 10*log10(sum(abs(x2_use).^2)/sum(abs(x2_use - x2_est_aligned).^2));

    x2_demod = qamdemod(x2_use,Q,'UnitAveragePower',true);

    x2_est_demod = qamdemod(x2_est_aligned,Q,'UnitAveragePower',true);

    errors2 = sum(x2_demod ~= x2_est_demod);

    %% Store results

    res(case_index).x1 = x1_use;
    res(case_index).x2 = x2_use;

    res(case_index).x1_est = x1_est;
    res(case_index).x2_est = x2_est;

    res(case_index).x1_est_aligned = x1_est_aligned;
    res(case_index).x2_est_aligned = x2_est_aligned;

    res(case_index).SINDR1 = SINDR1;
    res(case_index).SINDR2 = SINDR2;

    res(case_index).errors1 = errors1;
    res(case_index).errors2 = errors2;

    res(case_index).a0 = a0;
    res(case_index).a1 = a1;
    res(case_index).a2 = a2;
    res(case_index).a3 = a3;

    %% Print results

    fprintf('SINDR for signal 1: %.2f dB\n',SINDR1);
    fprintf('Symbol errors for signal 1: %d\n',errors1);

    fprintf('SINDR for signal 2: %.2f dB\n',SINDR2);
    fprintf('Symbol errors for signal 2: %d\n',errors2);

    if CFO_ON == 1 && CFO_COMP == 1
        fprintf('True CFO: %.6e\n',wT);
        fprintf('Estimated CFO: %.6e\n',CFO_est);
    end

end

%% Save case results

res1 = res(1);
res2 = res(2);
res3 = res(3);
res4 = res(4);

%% Nonlinearity spectrum test with A1 = 0

%fprintf('\nNonlinearity spectrum test with A1 = 0\n');

A1_image = 0;

%% Transmitter for nonlinearity image test

x1_up = upsample(x1,M1);
x2_up = upsample(x2,M2);

y1 = A1_image*conv(x1_up,H1);
y2 = A2*conv(x2_up,H2);

y1 = y1((length(H1)-1)/2+1:end-(length(H1)-1)/2);
y2 = y2((length(H2)-1)/2+1:end-(length(H2)-1)/2);

if length(y1) > length(y2)
    y2(end+1:length(y1)) = 0;
else
    y1(end+1:length(y2)) = 0;
end

N_signal = length(y1);
n = (0:N_signal-1).';

factor1 = exp(1i*omega1*n);
factor2 = exp(1i*omega2*n);

y1 = y1.*factor1;
y2 = y2.*factor2;

ytx = y1 + y2;

ytx_up = upsample(ytx,M);

ytx_filtered = conv(ytx_up,H_LP);
ytx_filtered = ytx_filtered((length(H_LP)-1)/2+1:end-(length(H_LP)-1)/2);

N_upsampled = length(ytx_filtered);
m = (0:N_upsampled-1).';

carrier_up = sqrt(2)*exp(1i*omega_c*m);

ytx_modulated = ytx_filtered.*carrier_up;

r = real(ytx_modulated);

a0 = 0.01;
a1 = 1;
a2 = -0.2/max(abs(r));
a3 = 0.15/(max(abs(r))^2);

yrx_nonlinear = a0 + a1*r + a2*r.^2 + a3*r.^3;

figure;
pwelch(r,[],[],[],fs_RF,'centered');
title('RF spectrum before nonlinearity, A1 = 0');

figure;
pwelch(yrx_nonlinear,[],[],[],fs_RF,'centered');
title('RF spectrum after nonlinearity, A1 = 0');

%% Scatter plots

Ns = 3000;

ideal_const = qammod(0:Q-1,Q,'UnitAveragePower',true);

figure;

subplot(2,2,1);
plot(real(res1.x2_est(1:min(Ns,length(res1.x2_est)))),imag(res1.x2_est(1:min(Ns,length(res1.x2_est)))),'.r','MarkerSize',12);
hold on;
plot(real(ideal_const),imag(ideal_const),'.b','MarkerSize',12);
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('Nonlinearity only, no compensation','FontSize',8);
set(gca,'FontSize',12);

subplot(2,2,2);
plot(real(res2.x2_est_aligned(1:min(Ns,length(res2.x2_est_aligned)))),imag(res2.x2_est_aligned(1:min(Ns,length(res2.x2_est_aligned)))),'.r','MarkerSize',12);
hold on;
plot(real(ideal_const),imag(ideal_const),'.b','MarkerSize',12);
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('Nonlinearity only, compensated','FontSize',8);
set(gca,'FontSize',12);

subplot(2,2,3);
plot(real(res3.x2_est(1:min(Ns,length(res3.x2_est)))),imag(res3.x2_est(1:min(Ns,length(res3.x2_est)))),'.r','MarkerSize',12);
hold on;
plot(real(ideal_const),imag(ideal_const),'.b','MarkerSize',12);
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('All errors, no compensation','FontSize',8);
set(gca,'FontSize',12);

subplot(2,2,4);
plot(real(res4.x2_est_aligned(1:min(Ns,length(res4.x2_est_aligned)))),imag(res4.x2_est_aligned(1:min(Ns,length(res4.x2_est_aligned)))),'.r','MarkerSize',12);
hold on;
plot(real(ideal_const),imag(ideal_const),'.b','MarkerSize',12);
grid on;
axis equal;
xlabel('In-phase');
ylabel('Quadrature');
title('All errors, compensated','FontSize',8);
set(gca,'FontSize',12);

%% Spectrum plots

figure;

subplot(3,2,1);
pwelch(res1.yrx_before_bpf,[],[],[],fs_RF,'centered');
title('Nonlinearity only: RF before filtering');

subplot(3,2,2);
pwelch(res1.yrx_after_nonlinearity,[],[],[],fs_RF,'centered');
title('Nonlinearity only: RF after nonlinearity');

subplot(3,2,3);
pwelch(res3.yrx_before_bpf,[],[],[],fs_RF,'centered');
title('All errors: before bandpass filter');

subplot(3,2,4);
pwelch(res4.yrx_after_bpf,[],[],[],fs_RF,'centered');
title('All errors: after bandpass filter');

subplot(3,2,5);
pwelch(res4.yrx_before_linearization,[],[],[],fs,'centered');
title('Baseband before linearization');

subplot(3,2,6);
pwelch(res4.yrx_after_linearization,[],[],[],fs,'centered');
title('Baseband after linearization');

%% Equalizer function

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