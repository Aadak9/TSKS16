%Matlab functions: conv, downsample, fftshift, find, freqz, qammod, qamdemod,
%rcosdesign, scatterplot, upsample.
B1 = 5e6;
B1_center = -5e6;
B2 = 10e6;
B2_center = -2.5e6;
fs = 30e6;
w1 = -fs*pi/3;
w2 = fs*pi/6;
A1 = 1;
A2 = 1;
M1 = 4;
M2 = 8;
L1 = 2^13;
L2 = 2^12;
S1 = 10;
S2 = 10;
rolloff1 = 1/3;
rolloff2 = 1/3;

g1=rcosdesign(rolloff1, S1, M1, 'sqrt')/sqrt(M1);
g2 = rcosdesign(rolloff1, S2, M2, 'sqrt')/sqrt(M2);

figure(1)

freqz(g1, L1); 
hold on;
freqz(g2, L2);
legend("g1", "g2");


figure(2)

subplot(2,2,1)
impz(g1.*g1)
title('Filter g1 cascade Impulse Response')

subplot(2,2,2)
impz(g2.*g2)
title('Filter g2 cascade Impulse Response')

subplot(2,2,3)
impz(downsample(M1*g1.*g1, M1))
title('G1 downsampled by M1')

subplot(2,2,4)
impz(downsample(M2*g2.*g2, M2))
title('G2 downsampled by M2')

%%
S1 = 1:40;
S2 = 1:40;

cases = [1 1;     
         1 0.1;   
         0.1 1];  

SIDR1_all = zeros(3,40);
SIDR2_all = zeros(3,40);
error_x1_all = zeros(3,40);
error_x2_all = zeros(3,40);

for c = 1:size(cases,1)
    A1 = cases(c,1);
    A2 = cases(c,2);

    for Q = [4, 64]

        x1 = qammod(randi([0 Q-1], L1, 1)', Q);
        x2 = qammod(randi([0 Q-1], L2, 1)', Q);

        SIDR1 = zeros(1, length(S1));
        SIDR2 = zeros(1, length(S2));

        error_x1 = zeros(1, length(S1));
        error_x2 = zeros(1, length(S2));

        for S = S1

            % Transmitter
            g1 = rcosdesign(rolloff1, S, M1, 'sqrt')/sqrt(M1);
            g2 = rcosdesign(rolloff1, S, M2, 'sqrt')/sqrt(M2);
            H1 = M1*g1;
            H2 = M2*g2;

            x1_transmit = upsample(x1, M1);
            x2_transmit = upsample(x2, M2);

            x1_transmit = conv(x1_transmit, H1);
            x2_transmit = conv(x2_transmit, H2);

            x1_transmit = x1_transmit .* A1;
            x2_transmit = x2_transmit .* A2;

            n1 = 1:length(x1_transmit);
            n2 = 1:length(x2_transmit);

            x1_transmit = x1_transmit .* exp(1j * n1 * w1 / fs);
            x2_transmit = x2_transmit .* exp(1j * n2 * w2 / fs);

            num_of_zeros = length(x2_transmit) - length(x1_transmit);
            x1_transmit = [x1_transmit zeros(1, num_of_zeros)];
            y_tx = x1_transmit + x2_transmit;

            % Receiver
            n = 1:length(y_tx);
            x1_received = y_tx .* exp(-1j * n * w1 / fs);
            x2_received = y_tx .* exp(-1j * n * w2 / fs);

            x1_received = conv(x1_received, g1);
            x2_received = conv(x2_received, g2);

            x1_received = downsample(x1_received, M1);
            x2_received = downsample(x2_received, M2);

            x1_received = x1_received .* 1/A1;
            x2_received = x2_received .* 1/A2;

            x1_est = x1_received(S+1:end-S*2);
            x2_est = x2_received(S+1:end-S);

            SIDR1(S) = 10*log10(sum(abs(x1).^2)/sum(abs(x1_est - x1).^2));
            SIDR2(S) = 10*log10(sum(abs(x2).^2)/sum(abs(x2_est - x2).^2));

            error_x1(S) = length(find(qamdemod(x1, Q) - qamdemod(x1_est, Q)));
            error_x2(S) = length(find(qamdemod(x2, Q) - qamdemod(x2_est, Q)));

        end

        % Store result per case 
        SIDR1_all(c,:) = SIDR1;
        SIDR2_all(c,:) = SIDR2;
        error_x1_all(c,:) = error_x1;
        error_x2_all(c,:) = error_x2;


        if (A1 == 1 && A2 == 1)
            figure(4);
        elseif (A1 == 1 && A2 == 0.1)
            figure(6);
        elseif (A1 == 0.1 && A2 == 1)
            figure(8);
        end

        if Q == 4
            p = 1;
        else
            p = 2;
        end

        subplot(2,1,p)
        plot(S1, SIDR1); hold on
        plot(S1, SIDR2)
        xlabel('S')
        ylabel('SIDR (dB)')
        title(['A1=',num2str(A1),' A2=',num2str(A2),' Q=',num2str(Q)])
        legend('x1','x2')

        if (A1 == 1 && A2 == 1)
            figure(5);
        elseif (A1 == 1 && A2 == 0.1)
            figure(7);
        elseif (A1 == 0.1 && A2 == 1)
            figure(9);
        end

        subplot(2,1,p)
        plot(S1, error_x1); hold on
        plot(S1, error_x2)
        xlabel('S')
        ylabel('Errors')
        title(['A1=',num2str(A1),' A2=',num2str(A2),' Q=',num2str(Q)])
        legend('x1','x2')

    end
end

%% --- A1 = 1, S1 = S2 = 20, A2 sweep 

A1 = 1;
S1 = 20;
S2 = 20;

A2_vec = logspace(0, -6, 120);

figure;
hold on;
grid on;
set(gca,'XScale','log');

for Q = [4, 64]

    x1 = qammod(randi([0 Q-1], L1, 1)', Q);
    x2 = qammod(randi([0 Q-1], L2, 1)', Q);

    errors_x1 = zeros(size(A2_vec));
    errors_x2 = zeros(size(A2_vec));

    for k = 1:length(A2_vec)

        A2 = A2_vec(k);

        
        g1 = rcosdesign(rolloff1, S1, M1, 'sqrt')/sqrt(M1);
        g2 = rcosdesign(rolloff2, S2, M2, 'sqrt')/sqrt(M2);

        H1 = M1*g1;
        H2 = M2*g2;

        %Transmitter 
        x1_transmit = upsample(x1, M1);
        x2_transmit = upsample(x2, M2);

        x1_transmit = conv(x1_transmit, H1);
        x2_transmit = conv(x2_transmit, H2);

        x1_transmit = x1_transmit .* A1;
        x2_transmit = x2_transmit .* A2;

        n1 = 1:length(x1_transmit);
        n2 = 1:length(x2_transmit);

        x1_transmit = x1_transmit .* exp(1j * n1 * w1 / fs);
        x2_transmit = x2_transmit .* exp(1j * n2 * w2 / fs);

        num_of_zeros = length(x2_transmit) - length(x1_transmit);

        if num_of_zeros > 0
            x1_transmit = [x1_transmit zeros(1, num_of_zeros)];
        elseif num_of_zeros < 0
            x2_transmit = [x2_transmit zeros(1, -num_of_zeros)];
        end

        y_tx = x1_transmit + x2_transmit;

        % Receiver x1
        n = 1:length(y_tx);

        x1_received = y_tx .* exp(-1j * n * w1 / fs);
        x1_received = conv(x1_received, g1);
        x1_received = downsample(x1_received, M1);
        x1_received = x1_received .* 1/A1;

        x1_est = x1_received(S1+1:end-2*S1);

        % Receiver x2
        x2_received = y_tx .* exp(-1j * n * w2 / fs);
        x2_received = conv(x2_received, g2);
        x2_received = downsample(x2_received, M2);
        x2_received = x2_received .* 1/A2;

        x2_est = x2_received(S2+1:end-S2);

        
        L = min([length(x1_est), length(x1), length(x2_est), length(x2)]);

        errors_x1(k) = sum(qamdemod(x1(1:L),Q) ~= qamdemod(x1_est(1:L),Q));
        errors_x2(k) = sum(qamdemod(x2(1:L),Q) ~= qamdemod(x2_est(1:L),Q));

    end

    %Threshold detection 
    idx = find(errors_x1 == 0 & errors_x2 == 0, 1, 'last');

    if isempty(idx)
        fprintf('Q = %d: No error-free A2 found\n', Q);
    else
        fprintf('Q = %d: Smallest A2 ≈ %.5g\n', Q, A2_vec(idx));
    end

    semilogx(A2_vec, errors_x1, '-o');
    semilogx(A2_vec, errors_x2, '-x');

end

xlabel('A2');
ylabel('Symbol Errors');
title('A1 = 1, S1 = S2 = 20');

legend('x1 (Q=4)','x2 (Q=4)','x1 (Q=64)','x2 (Q=64)');