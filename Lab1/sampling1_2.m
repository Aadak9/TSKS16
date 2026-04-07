clear; clc; close all;

%% =========================
% Parameters
% =========================
fs = 40e6;                 
T = 1/fs;
N = 8192;                 
n = 0:N-1;
t = n*T;

fc1 = 10e6;               % Desired center frequency
fc2 = 30.5e6;             % Undesired center frequency
Delta = 1/8;              

k_vals = -5:5;

A2 = 1;

P_vals = 1:40;            % Filter orders

A1_vals = [1, 0.01];
Amax_vals = [0.1, 0.01];  % dB

%% =========================
% Generate analog signals
% =========================

xa1 = zeros(1,N);
xa2 = zeros(1,N);

for k = k_vals
    f1k = fc1 + fc1 * k * Delta;
    f2k = fc2 + fc1 * k * Delta;

    xa1 = xa1 + sin(2*pi*f1k*t);
    xa2 = xa2 + sin(2*pi*f2k*t);
end

%% =========================
% Frequency axis for FFT
% =========================
w = 2*pi*(0:N-1)/N;          % digital frequency
omega = w * fs;              % analog frequency (rad/s)

wc = 0.8125 * pi * fs;       % cutoff frequency

%% =========================
% Main SIDR computation
% =========================

for A1 = A1_vals
for Amax = Amax_vals

    epsilon = sqrt(10^(Amax/10) - 1);

    SIDR = zeros(size(P_vals));

    for idx = 1:length(P_vals)
        P = P_vals(idx);

        % Full signal
        x = A1*xa1 + A2*xa2;

        % FFT
        X = fft(x);

        % Butterworth magnitude response
        H = 1 ./ sqrt(1 + epsilon^2 * (omega./wc).^(2*P));

        % Apply filter
        Xf = X .* H;

        % Back to time domain
        xf = real(ifft(Xf));

        % Desired signal (reference)
        x0 = A1 * xa1;

        % SIDR
        SIDR(idx) = 10*log10( sum(x0.^2) / sum((xf - x0).^2) );
    end

    %% Plot SIDR
    figure;
    plot(P_vals, SIDR, 'LineWidth', 1.5);
    grid on;
    xlabel('Filter order P');
    ylabel('SIDR (dB)');
    title(['SIDR vs P | A1 = ', num2str(A1), ...
           ', Amax = ', num2str(Amax), ' dB']);

end
end

%% =========================
% Example: Spectrum + Filter
% =========================

P = 20;              % Choose representative order
A1 = 1;
Amax = 0.1;

epsilon = sqrt(10^(Amax/10) - 1);

x = A1*xa1 + A2*xa2;

X = fft(x);

H = 1 ./ sqrt(1 + epsilon^2 * (omega./wc).^(2*P));

Xf = X .* H;
xf = real(ifft(Xf));

%% Apply Blackman-Harris window
wBH = blackmanharris(N)';
Xw = fft(xf .* wBH);
A_dB = 20*log10(abs(Xw)/max(abs(Xw)));

%% Frequency axis (normalized)
w_norm = w/pi;

%% Plot spectrum
figure;
plot(w_norm, A_dB, 'LineWidth', 1);
grid on;
xlim([0 1]);
ylim([-150 0]);
xlabel('Normalized frequency (\omegaT/\pi)');
ylabel('Magnitude (dB)');
title(['Spectrum after filtering (P = ', num2str(P), ')']);

%% Plot filter response
figure;
plot(w_norm, 20*log10(abs(H)/max(abs(H))), 'LineWidth', 1);
grid on;
xlabel('Normalized frequency (\omegaT/\pi)');
ylabel('Magnitude (dB)');
title(['Butterworth Filter Response (P = ', num2str(P), ')']);