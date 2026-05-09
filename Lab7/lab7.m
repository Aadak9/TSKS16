clc;
clear;
close all;

%% General simulation parameters

fs = 30e6;
fc = 350e6;

wc = 2*pi*fc;

SNR_dB = 30;

Mmod = 16;

N = 512;
Neq = 24;

%% Transmitter parameters

L1 = 2^13;
L2 = 2^12;

M = 100;

M1 = 4;
M2 = 8;

A1 = 0.1;
A2 = 1;

omega1 = -pi/3;
omega2 = pi/6;

rolloff = 1/3;

S1 = 16;
S2 = 16;
S = 10;

%% IQ mismatch parameters

g1 = 0.99;
g2 = 1.02;

phi1 = 0.5;
phi2 = 0.55;

%% CFO / PO parameters

w0T = 5*pi*1e-5;
alpha = 0.1*pi;

%% SFM parameter

r = 3e-5;

%% Channel

c = 0.25*[1 zeros(1,15) 2.4 zeros(1,15) 1];

%% Experiment control

IQM_ON = 1;
IQM_COMP = 1;

CFO_ON = 1;
CFO_COMP = 1;

SFM_ON = 1;

CHANNEL_ON = 1;
EQ_ON = 1;

%% Signal generation

x1 = qammod(randi([0 Mmod-1],L1,1).',Mmod);
x2 = qammod(randi([0 Mmod-1],L2,1).',Mmod);

%% Pulse shaping filters

G1 = rcosdesign(rolloff,S1,M1,"sqrt")/sqrt(M1);
G2 = rcosdesign(rolloff,S2,M2,"sqrt")/sqrt(M2);

H1 = M1*G1;
H2 = M2*G2;

%% Transmitter

x1_tx = upsample(x1,M1);
x1_tx = conv(x1_tx,H1);
x1_tx = A1*x1_tx;

x2_tx = upsample(x2,M2);
x2_tx = conv(x2_tx,H2);
x2_tx = A2*x2_tx;

n1 = 0:length(x1_tx)-1;
n2 = 0:length(x2_tx)-1;

x1_tx = x1_tx .* exp(1j*omega1*n1);
x2_tx = x2_tx .* exp(1j*omega2*n2);

L = max(length(x1_tx),length(x2_tx));

x1_tx = [x1_tx zeros(1,L-length(x1_tx))];
x2_tx = [x2_tx zeros(1,L-length(x2_tx))];

y_tx = x1_tx + x2_tx;

%% Add redundancy for CFO estimation

if CFO_ON
    y_tx = [y_tx(1:N) y_tx];
end

%% RF upconversion

G_lp = rcosdesign(rolloff,S,M,"sqrt")/sqrt(M);
H_lp = M*G_lp;

y_rx = upsample(y_tx,M);

y_rx = conv(y_rx,H_lp);

m = 0:length(y_rx)-1;

y_rx = y_rx .* (sqrt(2)*exp(1j*wc/(M*fs)*m));

y_rx = real(y_rx);

%% Channel

if CHANNEL_ON
    y_rx = conv(y_rx,c);
end

%% AWGN

y_rx = awgn(y_rx,SNR_dB,'measured');

%% Receiver downconversion

m = 0:length(y_rx)-1;

TH = 1/(M*fs);

if IQM_ON && CFO_ON

    lo = g1*cos((wc-w0T)*TH*m + phi1 - alpha) - 1j*g2*sin((wc-w0T)*TH*m + phi2 - alpha);

elseif IQM_ON && ~CFO_ON

    lo = g1*cos(wc*TH*m + phi1) - 1j*g2*sin(wc*TH*m + phi2);

elseif ~IQM_ON && CFO_ON

    lo = exp(-1j*((wc-w0T)*TH*m - alpha));

else

    lo = exp(-1j*wc*TH*m);

end

y_rx = y_rx .* (sqrt(2)*lo);

%% Store IQM distorted signal

y_iqm_before = y_rx;

%% Matched filter

y_rx = conv(y_rx,G_lp);

%% Downsample

y_rx = downsample(y_rx,M);

%% Remove filter delay

if CHANNEL_ON
    y_rx = y_rx(S+1:end-S-1);
else
    y_rx = y_rx(S+1:end-S);
end

%% Apply SFM

if SFM_ON
    y_rx = Interpolation_Farrow(y_rx,1+r);
end

%% IQM estimation and compensation (Algorithm 2)

if IQM_ON

    yr = real(y_rx);
    yi = imag(y_rx);

    Pr = mean(yr.^2);
    Pi = mean(yi.^2);

    g_est = sqrt(Pi/Pr);

    yp = -yr .* yi;
    Pyp = mean(yp);

    phi_est = Pyp / Pr;

    yr_comp = yr * g_est;
    yi_comp = yi + yr_comp * phi_est;

    if IQM_COMP
        y_rx = yr_comp + 1j*yi_comp;
    end

end

%% Store compensated IQ signal

y_iq_comp = y_rx;

%% CFO estimation

if CFO_ON

    w0T_est = angle(sum(conj(y_rx(1:N)) .* y_rx(N+1:2*N))) / N;

    fprintf('\nTrue CFO = %.6e\n',w0T);
    fprintf('Estimated CFO = %.6e\n',w0T_est);

end

%% CFO compensation

if CFO_COMP

    n = 1:length(y_rx);

    y_rx = y_rx .* exp(-1j*w0T_est*n);

end

%% Remove redundancy

if CFO_ON
    y_rx = y_rx(N+1:end);
end

%% Equalization

if EQ_ON

    if CFO_ON

        y_tx_ref = y_tx(N+1:end);

        if length(y_rx) > length(y_tx_ref)
            y_rx = y_rx(1:length(y_tx_ref));
        else
            y_tx_ref = y_tx_ref(1:length(y_rx));
        end

        y_end = length(y_tx_ref);

    else

        y_tx_ref = y_tx;

        if length(y_rx) > length(y_tx_ref)
            y_rx = y_rx(1:length(y_tx_ref));
        else
            y_tx_ref = y_tx_ref(1:length(y_rx));
        end

        y_end = length(y_tx_ref);

    end

    [h_eq,d_min,~,y_eq] = fir_eq(y_tx_ref,y_rx,Neq);

    y_comp = y_eq(d_min+1:y_end+d_min);

else

    y_comp = y_rx;

end

%% User 1 receiver

n = 0:length(y_comp)-1;

x1_rx = y_comp .* exp(-1j*omega1*n);

x1_rx = conv(x1_rx,G1);

x1_rx = downsample(x1_rx,M1);

x1_rx = x1_rx / A1;

x1_est = x1_rx(S1+1:S1+L1);

%% User 2 receiver

x2_rx = y_comp .* exp(-1j*omega2*n);

x2_rx = conv(x2_rx,G2);

x2_rx = downsample(x2_rx,M2);

x2_rx = x2_rx / A2;

x2_est = x2_rx(S2+1:S2+L2);

%% Match lengths

Lx1 = min(length(x1),length(x1_est));
Lx2 = min(length(x2),length(x2_est));

x1 = x1(1:Lx1);
x1_est = x1_est(1:Lx1);

x2 = x2(1:Lx2);
x2_est = x2_est(1:Lx2);

%% SINDR calculation

x1_SINDR = 10*log10(sum(abs(x1).^2) / sum(abs(x1 - x1_est).^2));

x2_SINDR = 10*log10(sum(abs(x2).^2) / sum(abs(x2 - x2_est).^2));

%% Symbol errors

x1_err = sum(qamdemod(x1,Mmod) ~= qamdemod(x1_est,Mmod));

x2_err = sum(qamdemod(x2,Mmod) ~= qamdemod(x2_est,Mmod));

%% Results

fprintf('\nUser 1\n');
fprintf('SINDR = %.2f dB\n',x1_SINDR);
fprintf('Symbol errors = %d\n',x1_err);

fprintf('\nUser 2\n');
fprintf('SINDR = %.2f dB\n',x2_SINDR);
fprintf('Symbol errors = %d\n',x2_err);

%% Scatter plots

Ns = 3000;

figure;

subplot(2,2,1)

scatter(real(y_iqm_before(1:Ns)),imag(y_iqm_before(1:Ns)),'.')

grid on
axis equal

title('Before IQ Compensation')

subplot(2,2,2)

scatter(real(y_iq_comp(1:Ns)),imag(y_iq_comp(1:Ns)),'.')

grid on
axis equal

title('After IQ Compensation')

subplot(2,2,3)

scatter(real(y_comp(1:Ns)),imag(y_comp(1:Ns)),'.')

grid on
axis equal

title('After Equalization')

subplot(2,2,4)

scatter(real(x2_est(1:Ns)),imag(x2_est(1:Ns)),'.')

hold on

scatter(real(x2(1:Ns)),imag(x2(1:Ns)),'r.')

grid on
axis equal

legend('Estimated symbols','Ideal symbols')

title('User 2 Constellation')


%% Spectrum plots

figure

subplot(3,1,1)

pwelch(y_iqm_before,[],[],[],fs,'centered')

title('Spectrum Before IQ Compensation')

subplot(3,1,2)

pwelch(y_iq_comp,[],[],[],fs,'centered')

title('Spectrum After IQ Compensation')

subplot(3,1,3)

pwelch(y_comp,[],[],[],fs,'centered')

title('Spectrum After Equalization')

%% SFM study

r_vec = [3e-3 3e-4 3e-5 3e-6 3e-7 3e-8];

SINDR_vec = zeros(size(r_vec));

for k = 1:length(r_vec)

    y_sfm = Interpolation_Farrow(y_tx,1+r_vec(k));

    Ls = min(length(y_sfm),length(y_tx));

    y_sfm = y_sfm(1:Ls);

    y_ref = y_tx(1:Ls);

    SINDR_vec(k) = 10*log10(sum(abs(y_ref).^2) / sum(abs(y_ref-y_sfm).^2));

end

figure

semilogx(r_vec,SINDR_vec,'-o')

grid on

xlabel('r')

ylabel('SINDR (dB)')

title('Effect of Sampling Frequency Mismatch')

%% Equalizer function

function [h_eq,d_min,error_min,y_rx_eq] = fir_eq(y_tx,y_rx,Neq)

d_min = 0;

error_min = Inf;

col = [y_rx(:); zeros(Neq,1)];
row = [y_rx(1) zeros(1,Neq)];

A = toeplitz(col,row);

B = (A'*A)\A';

for d = 0:Neq

    y_des = [zeros(1,d) y_tx zeros(1,Neq-d)];

    y_des = y_des(:);

    h_tmp = B*y_des;

    y_eq_tmp = A*h_tmp;

    err = mean(abs(y_eq_tmp-y_des).^2);

    if err < error_min

        error_min = err;

        h_eq = h_tmp;

        d_min = d;

        y_rx_eq = y_eq_tmp;

    end

end

y_rx_eq = y_rx_eq.';

end