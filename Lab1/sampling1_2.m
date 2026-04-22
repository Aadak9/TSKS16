clear; clc; close all;

fs = 40e6;
N = 8192;
t = (0:N-1)/fs;

fc1 = 10e6;
fc2 = 30.5e6;
Delta = 1/8;
k_vals = -5:5;

A2 = 1;
P_vals = 1:40;

A1_vals = [1, 0.01];
Amax_vals = [0.1, 0.01];
Delta_phase_vals = [0, 0.1, 0.01];

wc = 0.8125*pi*fs;

f = (0:N-1)*(fs/N);
omega = 2*pi*f;


%% 1. SIDR vs P


figure;

for d = 1:3

    phi = Delta_phase_vals(d);

    subplot(3,1,d);
    hold on; grid on;

    legends = {};

    for A1 = A1_vals
    for Amax = Amax_vals

        eps = sqrt(10^(Amax/10) - 1);
        SIDR = zeros(size(P_vals));

        for i = 1:length(P_vals)

            P = P_vals(i);

            xa1 = zeros(1,N);
            xa2 = zeros(1,N);
            xa1_ref = zeros(1,N);

            for k = k_vals

                f1 = fc1 + fc1*k*Delta;
                f2 = fc2 + fc1*k*Delta;

                ph = phi * (-1)^k;

                xa1 = xa1 + sin(2*pi*f1*t + ph);
                xa2 = xa2 + sin(2*pi*f2*t);

                xa1_ref = xa1_ref + sin(2*pi*f1*t);
            end

            x = A1*xa1 + A2*xa2;
            x0 = A1*xa1_ref;

            X = fft(x);

            
            H = 1 ./ sqrt(1 + eps^2 * ((omega)./wc).^(2*P));

            xf = real(ifft(X .* H));

            SIDR(i) = 10*log10(sum(x0.^2) / sum((xf - x0).^2));
        end

        plot(P_vals, SIDR, 'LineWidth', 1.5);

        legends{end+1} = sprintf('A1=%.2g, Amax=%.2g', A1, Amax);

    end
    end

    title(['SIDR vs P, phase = ', num2str(phi)]);
    xlabel('P');
    ylabel('SIDR (dB)');
    legend(legends,'Location','bestoutside');
end

%% =========================
% 2. Spectrum comparison
%% =========================

figure;

P_test = [10 20 40 60];

for d = 1:3

    phi = Delta_phase_vals(d);

    subplot(3,1,d);
    hold on; 
    grid on;

    legends = {};

    for p = 1:length(P_test)

        P = P_test(p);

        xa1 = zeros(1,N);
        xa2 = zeros(1,N);

        for k = k_vals

            f1 = fc1 + fc1*k*Delta;
            f2 = fc2 + fc1*k*Delta;

            ph = phi * (-1)^k;

            xa1 = xa1 + sin(2*pi*f1*t + ph);
            xa2 = xa2 + sin(2*pi*f2*t);
        end

        x = A1_vals(1)*xa1 + A2*xa2;

        X = fftshift(fft(x));

        eps = sqrt(10^(0.1/10) - 1);

        H = 1 ./ sqrt(1 + eps^2 * ((omega)./wc).^(2*P));

        xf = real(ifft(ifftshift(X .* H)));

        Xw = fftshift(fft(xf .* blackmanharris(N)'));

        A_dB = 20*log10(abs(Xw)/max(abs(Xw)));

        plot(f - fs/2, A_dB, 'LineWidth', 1.5);

        legends{p} = ['P = ', num2str(P)];
    end

    xlim([-fs/2 fs/2]);
    ylim([-220 0]);

    title(['Spectrum, phase = ', num2str(phi)]);
    xlabel('Hz');
    ylabel('dB');
    legend(legends);
end


%% 3. Filter response


figure; 
hold on; 
grid on;

f_plot = linspace(0, fs, 3000);
w_plot = 2*pi*f_plot;

eps = sqrt(10^(0.1/10) - 1);

for P = [10 20 40 60]

    H = 1 ./ sqrt(1 + eps^2 * (w_plot./wc).^(2*P));

    plot(f_plot, H./max(H), 'LineWidth', 1.5, ...
        'DisplayName', ['P = ', num2str(P)]);
end

ylim([0 1.05]);
xlabel('Hz');
ylabel('Normalized gain');
title('Filter response');
legend;


