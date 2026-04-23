clear; clc; close all;

N = 8192;
n = 0:N-1;
f = (0:N-1)/N;

w_axis = (0:N-1)*(2*pi/N);

A_list = [0.01 0.001];   % two amplitudes 

w_rect = ones(1,N);
w_bh = blackmanharris(N)';

w1_base = 0.25*pi;

for aIdx = 1:length(A_list)

    A2 = A_list(aIdx);

    figure;

    
    % CASE 1: RECT - COHERENT
    
    p = 2;
    dw = 2*pi*p/N;

    w1 = w1_base;
    w2 = w1 + dw;

    x = sin(w1*n) + A2*sin(w2*n);

    X = fft(x .* w_rect);
    X = 20*log10(abs(X)/max(abs(X)));

    subplot(2,2,1);
    plot(w_axis/pi, X);
    title(['Rect Coherent, A = ', num2str(A2)]);
    xlabel('\omega/\pi'); ylabel('dB'); grid on;

    % CASE 2: RECT - NONCOHERENT
    p = 700;
    dw = 2*pi*p/N;

    w1 = (0.25 + 1/N)*pi;
    w2 = w1 + dw;

    x = sin(w1*n) + A2*sin(w2*n);

    X = fft(x .* w_rect);
    X = 20*log10(abs(X)/max(abs(X)));

    subplot(2,2,2);
    plot(w_axis/pi, X);
    title(['Rect Noncoherent, A = ', num2str(A2)]);
    xlabel('\omega/\pi'); ylabel('dB'); grid on;


    % CASE 3: BH - COHERENT
    p = 5;
    dw = 2*pi*p/N;

    w1 = w1_base;
    w2 = w1 + dw;

    x = sin(w1*n) + A2*sin(w2*n);

    X = fft(x .* w_bh);
    X = 20*log10(abs(X)/max(abs(X)));

    subplot(2,2,3);
    plot(w_axis/pi, X);
    title(['BH Coherent, A = ', num2str(A2)]);
    xlabel('\omega/\pi'); ylabel('dB'); grid on;


    % CASE 4: BH - NONCOHERENT

    p = 5;
    dw = 2*pi*p/N;

    w1 = (0.25 + 1/N)*pi;
    w2 = w1 + dw;

    x = sin(w1*n) + A2*sin(w2*n);

    X = fft(x .* w_bh);
    X = 20*log10(abs(X)/max(abs(X)));

    subplot(2,2,4);
    plot(w_axis/pi, X);
    title(['BH Noncoherent, A = ', num2str(A2)]);
    xlabel('\omega/\pi'); ylabel('dB'); grid on;

end
