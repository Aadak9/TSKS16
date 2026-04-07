clear; clc; close all;

%% Parameters
fs = 40e6;
T = 1/fs;
N = 8192;
n = 0:N-1;
t = n*T;

fc1 = 10e6;
fc2 = 30.5e6;
Delta = 1/8;

k_vals = -5:5;

A2 = 1;
P_vals = 1:40;

A1_vals = [1, 0.01];
Amax_vals = [0.1, 0.01];

Delta_phase_vals = [0.1, 0.01];

%% Frequency axis
w = 2*pi*(0:N-1)/N;
omega = w * fs;
wc = 0.8125*pi*fs;

%% Loop over cases
for A1 = A1_vals
for Amax = Amax_vals
for Delta_phase = Delta_phase_vals

    epsilon = sqrt(10^(Amax/10) - 1);
    SIDR = zeros(size(P_vals));

    for idx = 1:length(P_vals)
        P = P_vals(idx);

        % Initialize signals
        xa1 = zeros(1,N);
        xa2 = zeros(1,N);

        % Build signals with phase distortion on xa1
        for k = k_vals
            f1k = fc1 + fc1*k*Delta;
            f2k = fc2 + fc1*k*Delta;

            % Phase distortion ONLY for xa1
            phi_k = Delta_phase * (-1)^k;

            xa1 = xa1 + sin(2*pi*f1k*t + phi_k);
            xa2 = xa2 + sin(2*pi*f2k*t); % no phase distortion
        end

        % Full signal
        x = A1*xa1 + A2*xa2;

        % FFT
        X = fft(x);

        % Butterworth magnitude response
        H = 1 ./ sqrt(1 + epsilon^2 * (omega./wc).^(2*P));

        % Apply magnitude-only filter
        Xf = X .* H;

        % Back to time domain
        xf = real(ifft(Xf));

        % Desired reference (NO phase distortion!)
        xa1_ref = zeros(1,N);
        for k = k_vals
            f1k = fc1 + fc1*k*Delta;
            xa1_ref = xa1_ref + sin(2*pi*f1k*t);
        end
        x0 = A1 * xa1_ref;

        % SIDR
        SIDR(idx) = 10*log10( sum(x0.^2) / sum((xf - x0).^2) );
    end

    %% Plot
    figure;
    plot(P_vals, SIDR, 'LineWidth', 1.5);
    grid on;
    xlabel('Filter order P');
    ylabel('SIDR (dB)');
    title(['SIDR vs P | A1=', num2str(A1), ...
           ', Amax=', num2str(Amax), ...
           ', Phase \Delta=', num2str(Delta_phase)]);

end
end
end