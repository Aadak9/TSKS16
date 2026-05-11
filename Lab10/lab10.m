clear; clc; close all;

%% PART 1: Spectrum estimation

N = 64;
P = 16*N;
S = [6:26 38:57];
Q = 64;
Navg = 1000;

PSD = zeros(1,P);
PSD_expected = zeros(1,P);

for i = 1:Navg

    X = zeros(1,N);
    X(S) = qammod(randi([0 Q-1],1,length(S)),Q,'UnitAveragePower',true);

    x = ifft(X,N);
    Xp = fft(x,P);

    PSD = PSD + abs(Xp).^2;

end

PSD = PSD/Navg;
PSD_dB = 10*log10(PSD/max(PSD) + eps);

n = 0:N-1;

for k = S

    hk = (1/N)*exp(1i*2*pi*(k-1)*n/N);
    Hk = fft(hk,P);

    PSD_expected = PSD_expected + abs(Hk).^2;

end

PSD_expected_dB = 10*log10(PSD_expected/max(PSD_expected) + eps);

f = (0:P-1)/P;

figure;
plot(f,PSD_dB); hold on;
plot(f,PSD_expected_dB,'--');
grid on;
xlabel('Normalized frequency');
ylabel('Power spectrum [dB]');
legend('Simulated','Expected');
title('OFDM spectrum');

fprintf("Part 1: %d OFDM symbols averaged\n\n",Navg);

%% PART 2: PAPR

N_values = [32 1024 8192];
Q_values = [4 64];
Navg = 100;

for Q = Q_values

    figure;

    for ni = 1:length(N_values)

        N = N_values(ni);

        PAPR_mean = zeros(1,N);
        PAPR_var = zeros(1,N);

        for K = 1:N

            PAPR = zeros(1,Navg);

            for i = 1:Navg

                X = zeros(1,N);
                X(1:K) = qammod(randi([0 Q-1],1,K),Q,'UnitAveragePower',true);

                x = ifft(X,N);

                PAPR(i) = max(abs(x).^2)/mean(abs(x).^2);

            end

            PAPR_mean(K) = mean(PAPR);
            PAPR_var(K) = var(PAPR);

        end

        K_values = 1:N;

        if Q == 4
            C = K_values;
        else
            C = 7*K_values/3;
        end

        subplot(3,2,2*ni-1);
        plot(K_values,10*log10(PAPR_mean)); hold on;
        plot(K_values,10*log10(C),'--');
        grid on;
        xlabel('K');
        ylabel('PAPR [dB]');
        legend('Simulated','Theoretical');
        title("PAPR, N = " + N + ", Q = " + Q);

        subplot(3,2,2*ni);
        plot(K_values,C./PAPR_mean);
        grid on;
        xlabel('K');
        ylabel('C / PAPR');
        title("Ratio, N = " + N + ", Q = " + Q);

        fprintf("Q = %d, N = %d, Navg = %d\n",Q,N,Navg);
        fprintf("PAPR at K=N = %.2f dB\n",10*log10(PAPR_mean(end)));
        fprintf("Variance at K=N = %.4f\n\n",PAPR_var(end));

    end
end