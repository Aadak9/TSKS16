clear; clc; close all;

N = 8192;
n = 0:N-1;

%{
p1 = 2;  % you can change this value

dw = 2*pi*p1/N;

% =========================
% CASE 1: Rectangular - Coherent
% =========================
w1 = 0.25*pi;
w2 = w1 + dw;

x = sin(w1*n) + 0.01*sin(w2*n);

w = ones(1,N);  % rectangular

xw = x .* w;

X = fft(xw);

X_mag = abs(X);
X_mag = X_mag / max(X_mag);
X_dB = 20*log10(X_mag);

%figure;
%plot(X_dB);
xlim([0 1500]);
ylim([-300 0]);
xlabel('Frequency bin k');
ylabel('Magnitude (dB)');
title('Rectangular Window - Coherent');
%}

% =========================
% CASE 2: Rectangular - Noncoherent
% =========================
%{
p2 = 700; %If lower the required +-0.5dB will be missed
dw = 2*pi*p2/N;
w1 = (0.25 + 1/N)*pi;
w2 = w1 + dw;

x = sin(w1*n) + 0.01*sin(w2*n);

w = ones(1,N);  % rectangular

xw = x .* w;

X = fft(xw);

X_mag = abs(X);
X_mag = X_mag / max(X_mag);
X_dB = 20*log10(X_mag);

figure;
plot(X_dB);
xlim([0 3000]);
ylim([-100 0]);
xlabel('Frequency bin k');
ylabel('Magnitude (dB)');
title('Rectangular Window - Noncoherent');
%}

%{
% =========================
% CASE 3: Blackman-Harris - Coherent
% =========================
p3 = 5;
dw = 2*pi*p3/N;
w1 = 0.25*pi;
w2 = w1 + dw;

x = sin(w1*n) + 0.01*sin(w2*n);

w = blackmanharris(N)';

xw = x .* w;

X = fft(xw);

X_mag = abs(X);
X_mag = X_mag / max(X_mag);
X_dB = 20*log10(X_mag);

figure;
plot(X_dB);
xlim([0 2000]);
ylim([-250 0]);
xlabel('Frequency bin k');
ylabel('Magnitude (dB)');
title('Blackman-Harris Window - Coherent');
%}

%
% =========================
% CASE 4: Blackman-Harris - Noncoherent
% =========================
p4 = 5;
dw = 2*pi*p4/N;
w1 = (0.25 + 1/N)*pi;
w2 = w1 + dw;

x = sin(w1*n) + 0.01*sin(w2*n);

w = blackmanharris(N)';

xw = x .* w;

X = fft(xw);

X_mag = abs(X);
X_mag = X_mag / max(X_mag);
X_dB = 20*log10(X_mag);

figure;
plot(X_dB);
xlim([0 1500]);
ylim([-300 0]);
xlabel('Frequency bin k');
ylabel('Magnitude (dB)');
title('Blackman-Harris Window - Noncoherent');
%}