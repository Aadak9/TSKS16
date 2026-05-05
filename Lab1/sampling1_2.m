clear; clc; close all;

fs = 40e6;
N = 8192;
t = (0:N-1)/fs;


fc1 = 10e6;
fc2 = 30.5e6;
Delta = 1/8;
k_vals = -5:5;
f = (0:N-1)*(fs/N);
omega = 2*pi*f;
omega_c = 0.8125*pi;   % digital rad/sample cutoff
w_axis  = linspace(-pi, pi, N);
A2 = 1;
P_vals = 0:60;

A1_vals = [1, 0.01];
Amax_vals = [0.1, 0.01];

Delta_phase_vals = [0, 0.1, 0.01];   

omega_c = 0.8125*pi*fs;


x_a1 = zeros(1,N);

for k = k_vals
    f1 = fc1 + fc1*k*Delta;
    x_a1 = x_a1 + sin(2*pi*f1*t);
end

figure;

for d = 1:3

    phi = Delta_phase_vals(d);

    subplot(3,1,d);
    hold on; 
    grid on;

    legends = {};   

    for A1 = A1_vals
    for Amax = Amax_vals

        eps = sqrt(10^(Amax/10) - 1);

        SIDR = zeros(size(P_vals));

        for i = 1:length(P_vals)

            P = P_vals(i);

            x_a1_filt = zeros(1,N);
            x_a2_filt = zeros(1,N);

            for k = k_vals

                f1 = fc1 + fc1*k*Delta;
                f2 = fc2 + fc1*k*Delta;

                ph = phi * (-1)^k;

                H1 = 1/sqrt(1 + eps^2*((2*pi*f1)/omega_c)^(2*P));
                H2 = 1/sqrt(1 + eps^2*((2*pi*f2)/omega_c)^(2*P));

                x_a1_filt = x_a1_filt + H1*sin(2*pi*f1*t + ph);
                x_a2_filt = x_a2_filt + H2*sin(2*pi*f2*t);

            end

            x  = A1*x_a1_filt + A2*x_a2_filt;
            x0 = A1*x_a1;

            signal_power = sum(x0.^2);
            error_power  = sum((x - x0).^2) + eps;

            SIDR(i) = 10*log10(signal_power / error_power);

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


%% 2. Spectrum comparison

figure;

window = blackmanharris(N)';
F = (-N/2:N/2-1)*(fs/N);

P_test = [10 20 40 60];

omega_c = 0.8125*pi*fs;
eps = sqrt(10^(0.1/10) - 1);

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

        X = fftshift(fft(x .* window));

        w = 2*pi*F;  

        
        H = 1 ./ sqrt(1 + eps^2 * (abs(w)/omega_c).^(2*P));

        Xf = X .* H;

        A_dB = 20*log10(abs(Xf)/max(abs(Xf)) + 1e-12);

        plot(F, A_dB, 'LineWidth', 1.5);

        legends{end+1} = ['P = ', num2str(P)];

    end

    %xlim([-fs/2 fs/2]);
    %ylim([-220 10]);

    title(['Spectrum, phase = ', num2str(phi)]);
    xlabel('Frequency (Hz)');
    ylabel('Magnitude (dB)');
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

    H = 1 ./ sqrt(1 + eps^2 * (abs(w_plot)./omega_c).^(2*P));

    plot(f_plot, H./max(H), 'LineWidth', 1.5, ...
        'DisplayName', ['P = ', num2str(P)]);
end

ylim([0 1.05]);
xlabel('Hz');
ylabel('Normalized gain');
title('Filter response');
legend;


