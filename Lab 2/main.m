clc; clear; close all;


%% PARAMETERS

Bvals = [1 2 4 8 16];


%% PART 1: HISTOGRAMS OF QUANTIZATION ERROR

N = 2^20;

p = 7;%A prime is choosen so the error is not periodic
n = 0:N-1;
w1T = 2*pi*p/N;

signals1 = {
    randn(1,N)
    rand(1,N) - 0.5
    cos(w1T*n)
    cos(0.25*pi*n)
};

titles1 = {
    'Gaussian noise'
    'Uniform noise'
    'Cosine (prime frequency)'
    'Cosine (0.25\pi)'
};

for s = 1:length(signals1)

    x = AGCunity(signals1{s});

    figure('Name', titles1{s});

    for i = 1:length(Bvals)

        xq = quant(x, Bvals(i));
        e = xq - x;

        subplot(2,3,i)
        histogram(e, 100, 'Normalization', 'pdf')
        title(['B = ' num2str(Bvals(i))])
        xlabel('Error'); ylabel('PDF')

    end
end


%% PART 2: SNR VS N
B = 16;
Nvals = round(logspace(2,5,50));
SNR_theory = 6.02*B + 1.76;

figure; hold on;
colors = lines(3);

for s = 1:3

    SNR_vals = zeros(size(Nvals));

    for k = 1:length(Nvals)

        N = Nvals(k);
        n = 0:N-1;

        % Select signal
        if s == 1
            x = cos(0.5*sqrt(3)*pi*n);
            name = '0.5\sqrt{3}\pi';
        elseif s == 2
            x = cos(0.25*pi*n);
            name = '0.25\pi';
        else
            x = cos(0.5*pi*n);
            name = '0.5\pi';
        end

        % Quantization
        xq = quant(x, B);
        e = xq - x;

        % SNR
        SNR_vals(k) = 10*log10(sum(x.^2) / sum(e.^2));

    end

    semilogx(Nvals, SNR_vals, 'LineWidth', 1.5, 'Color', colors(s,:));

end

% Theory line
yline(SNR_theory, '--k', 'LineWidth', 1.2);

grid on;
xlabel('N');
ylabel('SNR (dB)');
title('SNR vs N for Different Cosine Signals');
legend('0.5\sqrt{3}\pi','0.25\pi','0.5\pi','Theory','Location','best');
xlim([0 1e4]);
%% PART 3: SPECTRUM ANALYSIS WINDOWING
N = 8192;
B = 16;

n = 0:N-1;
f = (0:N-1)/N;

w_coherent    = 0.25*pi;
w_noncoherent = (0.25 + 1/N)*pi;

colors = lines(4);

figure; hold on;


% Rectangular window

w = rectwin(N)';

% Coherent
x = cos(w_coherent*n);
xq = quant(x, B);
X = fft(xq .* w, N);
mag = 20*log10(2*abs(X)/N + eps);
mag = mag - max(mag);
plot(f(1:N/2), mag(1:N/2), 'Color', colors(1,:), 'LineWidth', 1.2);

% Noncoherent
x = cos(w_noncoherent*n);
xq = quant(x, B);
X = fft(xq .* w, N);
mag = 20*log10(2*abs(X)/N + eps);
mag = mag - max(mag);
plot(f(1:N/2), mag(1:N/2), '--', 'Color', colors(2,:), 'LineWidth', 1.2);


% Blackman-Harris 

w = blackmanharris(N)';

% Coherent
x = cos(w_coherent*n);
xq = quant(x, B);
X = fft(xq .* w, N);
mag = 20*log10(2*abs(X)/N + eps);
mag = mag - max(mag);
plot(f(1:N/2), mag(1:N/2), 'Color', colors(3,:), 'LineWidth', 1.2);

% Noncoherent
x = cos(w_noncoherent*n);
xq = quant(x, B);
X = fft(xq .* w, N);
mag = 20*log10(2*abs(X)/N + eps);
mag = mag - max(mag);
plot(f(1:N/2), mag(1:N/2), '--', 'Color', colors(4,:), 'LineWidth', 1.2);


grid on;
xlabel('Normalized Frequency');
ylabel('Magnitude (dB)');
title('Spectrum Comparison: Coherent vs Noncoherent, Window Effects');


legend( ...
    'Rectangular - Coherent', ...
    'Rectangular - Noncoherent', ...
    'Blackman-Harris - Coherent', ...
    'Blackman-Harris - Noncoherent' ...
);


%% PART 4: FFT LENGTH COMPARISON
B = 16;
N1 = 8192;
N2 = 16*N1;

Ns = [N1 N2];

figure; hold on;

colors = lines(2);

for i = 1:2

    N = Ns(i);
    n = 0:N-1;

    x = cos(0.5*sqrt(3)*pi*n);
    xq = quant(x, B);

    w = blackmanharris(N)';
    W = sum(w);

    xw = xq .* w;

    X = fft(xw, N);

    % Proper scaling + normalization
    mag = 20*log10(2*abs(X)/W + eps);
    mag = mag - max(mag);

    % Frequency axis (0 → 1)
    f = (0:N-1)/N;

    plot(f, mag, 'LineWidth', 1.2, 'Color', colors(i,:));

end

grid on;
xlabel('Normalized Frequency (0 to 1)');
ylabel('Magnitude (dB)');
title('Spectrum Comparison: N = 8192 vs N = 131072');
legend('N = 8192', 'N = 1131072');
