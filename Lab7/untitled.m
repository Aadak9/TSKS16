close all;
clear;

% ----- CODE -----

% Design variables
M = 100;
M1 = 4;
M2 = 8;
L1 = 2^13;
L2 = 2^12;
fs = 30000000;
fc = 350000000;
rolloff = 1/3;
%A1 = 1;
%A2 = 1;
S_max = 40;
%S_LP = 20;
%Q = 64;
%noise_SNR = 20;
c = zeros(33,1);
c(1) = 0.25;
c(17) = 2.4 * 0.25;
c(33) = 0.25;
wT = 5e-7 * pi;
alpha = 0.1 * pi;
N_offset = 512;

% Variables for lab 7
g1 = 0.99;
g2 = 1.02;
phi1 = 0.5;
phi2 = 0.55;
A1 = 0.1;
A2 = 1;
noise_SNR = 30;
S1 = 16;
S2 = 16;
S_LP = 10;
Q = 16;
r = 0;

% ----- Task 1 -----

task1 = 1;
if task1 == 1
    fprintf('Task 1:\n');
    transmultiplexer(A1, A2, M1, M2, L1, L2, S1, S2, Q, rolloff, M, fs, fc, S_LP, c, noise_SNR, wT, alpha, N_offset, g1, g2, phi1, phi2, r, 0)
    transmultiplexer(A1, A2, M1, M2, L1, L2, S1, S2, Q, rolloff, M, fs, fc, S_LP, c, noise_SNR, wT, alpha, N_offset, g1, g2, phi1, phi2, r, 1)
    fprintf('\n');
end

% ----- Task 2 -----

task2 = 1;
if task2 == 1
    fprintf('Task 2:\n');
    r = 3e-5;
    %transmultiplexer(A1, A2, M1, M2, L1, L2, S1, S2, Q, rolloff, M, fs, fc, S_LP, c, noise_SNR, wT, alpha, N_offset, g1, g2, phi1, phi2, r, 0)
    transmultiplexer(A1, A2, M1, M2, L1, L2, S1, S2, Q, rolloff, M, fs, fc, S_LP, c, noise_SNR, wT, alpha, N_offset, g1, g2, phi1, phi2, r, 1)
    
    r = 1e-5;
    transmultiplexer(A1, A2, M1, M2, L1, L2, S1, S2, Q, rolloff, M, fs, fc, S_LP, c, noise_SNR, wT, alpha, N_offset, g1, g2, phi1, phi2, r, 1)
    
    r = 3e-6;
    transmultiplexer(A1, A2, M1, M2, L1, L2, S1, S2, Q, rolloff, M, fs, fc, S_LP, c, noise_SNR, wT, alpha, N_offset, g1, g2, phi1, phi2, r, 1)
end

% ----- FUNCTIONS -----

% Equalizer function
function yrx_eq = eq_filter(y_tx, y_rx, N_eq)

    % Ensure column vectors
    y_tx = y_tx(:);
    y_rx = y_rx(:);
    
    % Vectors to store errors
    ERR = zeros(N_eq + 1, 1);
    best_error = Inf;
    
    % Go through all possible delays
    for d = 0:N_eq

        % Create the toeplitz matrix
        A = toeplitz([y_rx; zeros(N_eq,1)], [y_rx(1); zeros(N_eq,1)]); 

        % Align y_tx with delay d
        y_tx_d = [zeros(d,1); y_tx; zeros(N_eq-d,1)];

        % Least squares solution
        h = (A' * A) \ (A' * y_tx_d);
        current_yrx_eq = A * h;

        % Compute squared error
        ERR(d + 1) = sum(abs(current_yrx_eq - y_tx_d).^2);

        % Update the error
        if ERR(d + 1) < best_error
            best_error = ERR(d + 1);
            h_eq = h;
            d_min = d;
        end
    end

    yrx_eq = conv(y_rx, h_eq);
    yrx_eq = yrx_eq(d_min+1 : end - (N_eq-d_min));
end

% Function for the transmultiplexer
function transmultiplexer(A1, A2, M1, M2, L1, L2, S1, S2, Q, rolloff, M, fs, fc, S_LP, c, noise_SNR, wT, alpha, N_offset, g1, g2, phi1, phi2, r, IQM_comp)

    % ----- CREATE SIGNALS -----

    % Generate two independent random QAM symbol sequences
    x1 = qammod(randi([0 Q-1], L1, 1), Q, 'UnitAveragePower', true);
    x2 = qammod(randi([0 Q-1], L2, 1), Q, 'UnitAveragePower', true);

    % ----- TRANSMITTER -----

    % Upsample the symbol sequences
    x1_up = upsample(x1, M1);
    x2_up = upsample(x2, M2);

    % Define the modulation frequencies
    omega1 = -pi/3;
    omega2 = pi/6;

    % Arrays for the upcomming SIDR values
    %SIDR1 = zeros(S_max, 1);
    %SIDR2 = zeros(S_max, 1);

    % Loop through different filterlengths
    %for S = 1:S_max

    % Design square-root raised cosine filters for both channels
    G1 = rcosdesign(rolloff, S1, M1, 'sqrt')/sqrt(M1);
    G2 = rcosdesign(rolloff, S2, M2, 'sqrt')/sqrt(M2);

    % Create transmit filters
    H1 = M1 .* G1;
    H2 = M2 .* G2;

    % Filter the upsampled signals
    y1 = A1 * conv(x1_up, H1);
    y2 = A2 * conv(x2_up, H2);
    y1 = y1((length(H1)-1)/2 + 1 : end - (length(H1)-1)/2);
    y2 = y2((length(H2)-1)/2 + 1 : end - (length(H2)-1)/2);

    % Pad y2 with zeros so that y1 and y2 have the same length
    if length(y1) > length(y2)
        y2(end+1:length(y1)) = 0;
    else
        y1(end+1:length(y2)) = 0;
    end

    % Modulate both signals onto different frequencies
    N = length(y1);
    n = (0:N-1).';
    factor1 = exp(1i * n * omega1);
    factor2 = exp(1i * n * omega2);
    y1 = y1 .* factor1;
    y2 = y2 .* factor2;

    % Combine the two signals into a single transmission signal
    ytx = y1 + y2;

    % Add reduntant samples for the CFO and PO
    ytx = [ytx(1:N_offset); ytx];

    % ----- ZERO-IF ARCHITECTURE -----

    % Create filters
    G_LP = rcosdesign(rolloff, S_LP, M, 'sqrt')/sqrt(M);
    H_LP = M .* G_LP;

    % Upsampling
    ytx_up = upsample(ytx, M);

    % Filter with H
    ytx_filtered = conv(ytx_up, H_LP);
    ytx_filtered = ytx_filtered((length(H_LP)-1)/2 + 1 : end - (length(H_LP)-1)/2);

    % Modulate up with the carrier frequency
    N_upsampled = length(ytx_filtered);
    m = (0:N_upsampled-1).';
    omega_c = 2*pi*fc/(M*fs); 
    carrier_up = sqrt(2) * exp(1i * omega_c * m);
    ytx_modulated = ytx_filtered .* carrier_up;

    % Take the real part of the signal
    ytx_real = real(ytx_modulated);

    % Add the channel
    y_after_channel = conv(ytx_real, c);
    y_after_channel = y_after_channel((length(c)-1)/2 + 1 : end - (length(c)-1)/2);

    if noise_SNR ~= 0

        % Add noise
        y_after_channel = awgn(y_after_channel, noise_SNR, 'measured');

    end

    % ----- INTRODUCTION OF CFO AND PO -----

    % Add CFO and PO
    %carrier_down = sqrt(2) * exp(-1i * ((omega_c - wT) * m - alpha));
    carrier_down = sqrt(2) * (g1 * cos(2*pi*fc/(M*fs)*m + phi1) - 1i * g2 * sin(2*pi*fc/(M*fs)*m + phi2));
    yrx_modulated = y_after_channel .* carrier_down;

    % ----- ZERO-IF RECEIVER -----
    
    % Filter with G
    yrx_filtered = conv(yrx_modulated, G_LP);
    yrx_filtered = yrx_filtered((length(G_LP)-1)/2 + 1 : end - (length(G_LP)-1)/2);

    % Downsample
    yrx = downsample(yrx_filtered, M);

    % ----- COMPENSATION FOR CFO AND IQM -----

    if IQM_comp == 1
        % Compensate for IQM
        yrx = IQM_compensation(yrx);
    end

    % Estimate CFO
    CFO_est = angle(sum(conj(yrx(1:N_offset)) .* yrx(N_offset+1:2*N_offset))) / N_offset;

    % Compensate for the CFO
    n_temp = (0:length(yrx)-1).';
    yrx_CFO = yrx .* exp(-1i * CFO_est * n_temp);
    
    % Remove the redundant samples
    yrx = yrx_CFO(N_offset+1:end);
    ytx = ytx(N_offset+1:end);

    % Apply SFM interpolation
    if r ~= 0
        yrx = Interpolation_Farrow(yrx.', 1 + r);
        yrx = yrx.';
    end

    % ----- EQUALIZER -----

    % Define filter order
    N_eq = 10;

    % Run the equalizer
    yrx_eq = eq_filter(ytx, yrx, N_eq);

    % ----- RECEIVER -----

    % Demodulate received signal using conjugate carriers
    y_rx1 = yrx_eq .* conj(factor1);
    y_rx2 = yrx_eq .* conj(factor2);

    % Apply matched filtering
    r1 = conv(y_rx1, G1);
    r2 = conv(y_rx2, G2);
    r1 = r1((length(G1)-1)/2 + 1 : end - (length(G1)-1)/2);
    r2 = r2((length(G2)-1)/2 + 1 : end - (length(G2)-1)/2);

    % Downsample the signals
    x1_est = downsample(r1, M1);
    x2_est = downsample(r2, M2);

    % Compensate for the amplitude
    x1_est = x1_est / A1;
    x2_est = x2_est / A2;

    % Calculate the current SIDR
    %SIDR1(S) = 10 * log10(sum(abs(x1).^2) / sum(abs(x1 - x1_est).^2));
    %SIDR2(S) = 10 * log10(sum(abs(x2).^2) / sum(abs(x2 - x2_est).^2));
    SINDR1 = 10 * log10(sum(abs(x1).^2) / sum(abs(x1 - x1_est).^2));
    SINDR2 = 10 * log10(sum(abs(x2).^2) / sum(abs(x2 - x2_est).^2));
    fprintf('SINDR for signal 1: %d\n', SINDR1);
    fprintf('SINDR for signal 2: %d\n', SINDR2);

    % Create scatterplots for a chosen amount of S values
    %if S == 10
    figure;
    sgtitle(sprintf('S1 = %d and S2 = %d', S1, S2));

    subplot(2,2,1);
    plot(real(x1), imag(x1), '.');
    ylim([-1.1 1.1]);
    xlim([-1.1 1.1]);
    title('$x_1$', 'Interpreter', 'latex');

    subplot(2,2,2);
    plot(real(x1_est), imag(x1_est), '.');
    ylim([-1.1 1.1]);
    xlim([-1.1 1.1]);
    title('$\hat{x}_1$', 'Interpreter', 'latex');

    subplot(2,2,3);
    plot(real(x2), imag(x2), '.');
    ylim([-1.1 1.1]);
    xlim([-1.1 1.1]);
    title('$x_2$', 'Interpreter', 'latex');

    subplot(2,2,4);
    plot(real(x2_est), imag(x2_est), '.');
    ylim([-1.1 1.1]);
    xlim([-1.1 1.1]);
    title('$\hat{x}_2$', 'Interpreter', 'latex');
    %end

    % Calculate errors
    %if S == 10
    % Demodulate the signals
    x1_demod = qamdemod(x1, Q, 'UnitAveragePower', true);
    x2_demod = qamdemod(x2, Q, 'UnitAveragePower', true);
    x1_est_demod = qamdemod(x1_est, Q, 'UnitAveragePower', true);
    x2_est_demod = qamdemod(x2_est, Q, 'UnitAveragePower', true);

    % Find all errors
    errors1 = find(x1_demod ~= x1_est_demod);
    errors2 = find(x2_demod ~= x2_est_demod);

    fprintf('Errors for x1: %d\n', max(errors1));
    fprintf('Errors for x2: %d\n', max(errors2));
    %end
    %end

    % Plot SIDR
    %figure;
    %plot(1:S_max, SIDR1, 1:S_max, SIDR2);
    %legend('SIDR1', 'SIDR2');
    %xlabel('Filterlength S');
    %ylabel('SIDR (dB)');
    %title(sprintf('SIDR vs S (Q = %d)', Q));
    %grid on;
end

% Function for IQM compensation
function y = IQM_compensation(x)

    % Estimation
    yr=real(x);
    yi=imag(x);
    Pr=mean(yr.^2);
    Pi=mean(yi.^2);
    g_est=sqrt(Pi/Pr);
    yp=-yr.*yi;
    Pyp=mean(yp);
    phi_est=Pyp/Pr;
    
    % Compensation
    yr=yr*g_est;
    yi=yi+yr*phi_est;
    y=yr+1j*yi;

end