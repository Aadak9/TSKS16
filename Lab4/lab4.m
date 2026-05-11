clear; close all; clc;

L0 = 2^12;
Neq = 20;

c = 0.25*exp(1i*0.1*pi)*[1 2.4 1];

%% Case 1 white Gaussian input

y_tx_1 = randn(1,L0);
y_tx_1 = y_tx_1 / sqrt(mean(abs(y_tx_1).^2));

y_rx_1_clean = conv(y_tx_1,c);
y_rx_1 = awgn(y_rx_1_clean,30,"measured");

[h_eq_1,d_min_1,error_min_1,Lerr_1] = fir_eq(y_tx_1,y_rx_1,Neq);

SNDR_1 = 10*log10(sum(abs(y_tx_1(1+d_min_1:Lerr_1+d_min_1)).^2) / error_min_1);

fprintf("Case 1 white input\n");
fprintf("d_min = %d\n",d_min_1);
fprintf("ERR_min = %.4e\n",error_min_1);
fprintf("SNDR = %.2f dB\n\n",SNDR_1);

%% Case 2 filtered white Gaussian input

b = fir1(30,0.5);

y_tx_2 = randn(1,L0);
y_tx_2 = filter(b,1,y_tx_2);
y_tx_2 = y_tx_2 / sqrt(mean(abs(y_tx_2).^2));

y_rx_2_clean = conv(y_tx_2,c);
y_rx_2 = awgn(y_rx_2_clean,30,"measured");

[h_eq_2,d_min_2,error_min_2,Lerr_2] = fir_eq(y_tx_2,y_rx_2,Neq);

SNDR_2 = 10*log10(sum(abs(y_tx_2(1+d_min_2:Lerr_2+d_min_2)).^2) / error_min_2);

fprintf("Case 2 filtered input\n");
fprintf("d_min = %d\n",d_min_2);
fprintf("ERR_min = %.4e\n",error_min_2);
fprintf("SNDR = %.2f dB\n\n",SNDR_2);

%% Magnitude responses

figure;
[Hc,w] = freqz(c,1,2048);
[H1,~] = freqz(h_eq_1,1,2048);
[H2,~] = freqz(h_eq_2,1,2048);

plot(w/pi,20*log10(abs(Hc)+eps),'LineWidth',1.3); hold on;
plot(w/pi,20*log10(abs(H1)+eps),'LineWidth',1.3);
plot(w/pi,20*log10(abs(H2)+eps),'LineWidth',1.3);
grid on;
xlabel("Normalized frequency \omega/\pi");
ylabel("Magnitude [dB]");
legend("Channel","Equalizer white","Equalizer filtered","Location","best");
title("Magnitude responses");

%% SNDR versus equalizer order

Neq_range = 1:50;

SNDR_white = zeros(size(Neq_range));
SNDR_filt = zeros(size(Neq_range));

d_white = zeros(size(Neq_range));
d_filt = zeros(size(Neq_range));

for k = 1:length(Neq_range)

    Neq = Neq_range(k);

    [~,d_white(k),ERR_white,Lerr_white] = fir_eq(y_tx_1,y_rx_1,Neq);
    [~,d_filt(k),ERR_filt,Lerr_filt] = fir_eq(y_tx_2,y_rx_2,Neq);

    P_white = sum(abs(y_tx_1(1+d_white(k):Lerr_white+d_white(k))).^2);
    P_filt = sum(abs(y_tx_2(1+d_filt(k):Lerr_filt+d_filt(k))).^2);

    SNDR_white(k) = 10*log10(P_white / ERR_white);
    SNDR_filt(k) = 10*log10(P_filt / ERR_filt);

end

figure;
plot(Neq_range,SNDR_white,'LineWidth',1.5); hold on;
plot(Neq_range,SNDR_filt,'LineWidth',1.5);
grid on;
xlabel("N_{eq}");
ylabel("SNDR [dB]");
legend("White input","Filtered input","Location","best");
title("SNDR versus equalizer order");

figure;
plot(Neq_range,d_white,'LineWidth',1.5); hold on;
plot(Neq_range,d_filt,'LineWidth',1.5);
grid on;
xlabel("N_{eq}");
ylabel("d_{min}");
legend("White input","Filtered input","Location","best");
title("Best delay versus equalizer order");

%% Equalizer function

function [h_eq,d_min,ERR_min,L_best] = fir_eq(y_tx,y_rx,Neq)

    % Convert signals to column vectors
    y_tx = y_tx(:);
    y_rx = y_rx(:);

    % Length of transmitted and received signals
    Ltx = length(y_tx);
    Lrx = length(y_rx);

    % Initialize minimum error
    ERR_min = Inf;

    % Initialize best delay
    d_min = 0;

    % Initialize equalizer coefficients
    h_eq = zeros(Neq+1,1);

    % Initialize best signal length
    L_best = 0;

    % Test all delays from 0 to Neq
    for d = 0:Neq

        % Number of valid samples
        L = min(Lrx-Neq,Ltx-d);

        % Build Toeplitz matrix from received signal
        %
        % Each row contains delayed samples:
        %
        % [y_rx(n) y_rx(n-1) ... y_rx(n-Neq)]
        %
        col = y_rx(Neq+1:Neq+L);
        row = y_rx(Neq+1:-1:1).';

        % Toeplitz convolution matrix
        A = toeplitz(col,row);

        % Desired delayed transmitted signal
        %
        % Equalizer output should approximate:
        %
        % y_tx(n-d)
        %
        b = y_tx(1+d:L+d);

        % Least-squares solution
        %
        % Minimize:
        %
        % ||A*h_eq - b||^2
        %
        h_tmp = (A'*A)\(A'*b);

        % Equalized received signal
        y_eq = A*h_tmp;

        % Compute squared error
        ERR = sum(abs(y_eq - b).^2);

        % Save result if error is smaller
        if ERR < ERR_min

            % Update minimum error
            ERR_min = ERR;

            % Save equalizer coefficients
            h_eq = h_tmp;

            % Save optimal delay
            d_min = d;

            % Save valid signal length
            L_best = L;
        end
    end
end