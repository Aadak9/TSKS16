clear; close all;

rolloff1 = 1/3;
rolloff2 = 1/3;

M1 = 4;
M2 = 8;

L1 = 2^13;
L2 = 2^12;

fs = 30e6;
omega1 = -pi/3;
omega2 = pi/6;

S1 = 1:40;

G1 = rcosdesign(rolloff1, 10, M1, 'sqrt')/sqrt(M1);
G2 = rcosdesign(rolloff2, 10, M2, 'sqrt')/sqrt(M2);

figure;
freqz(G1); hold on;
freqz(G2);
legend("G1","G2");

figure;
subplot(2,2,1); impz(G1.*G1); title("G1 cascaded");
subplot(2,2,2); impz(G2.*G2); title("G2 cascaded");
subplot(2,2,3); impz(downsample(M1*G1.*G1,M1)); title("M1G1G1 ↓M1");
subplot(2,2,4); impz(downsample(M2*G2.*G2,M2)); title("M2G2G2 ↓M2");

cases = [1 1;
         1 0.1;
         0.1 1];

SIDR1_all = zeros(3,40);
SIDR2_all = zeros(3,40);

y_store = [];

scatter_saved = false;

for c = 1:3

    A1 = cases(c,1);
    A2 = cases(c,2);

    for Q = [4, 64]

        x1 = qammod(randi([0 Q-1], L1, 1)', Q);
        x2 = qammod(randi([0 Q-1], L2, 1)', Q);

        SIDR1 = zeros(1,40);
        SIDR2 = zeros(1,40);
        err1 = zeros(1,40);
        err2 = zeros(1,40);

        for S = S1

            G1 = rcosdesign(rolloff1, S, M1, 'sqrt')/sqrt(M1);
            G2 = rcosdesign(rolloff2, S, M2, 'sqrt')/sqrt(M2);

            H1 = M1*G1;
            H2 = M2*G2;

            x1_tx = upsample(x1, M1);
            x1_tx = conv(x1_tx, H1) .* A1;
            n1_up = 1:length(x1_tx);
            x1_tx = x1_tx .* exp(1i*omega1*n1_up);

            x2_tx = upsample(x2, M2);
            x2_tx = conv(x2_tx, H2) * A2;
            x2_tx = x2_tx .* exp(1i*omega2*(1:length(x2_tx)));

            diff = length(x2_tx) - length(x1_tx);
            x1_tx = [x1_tx zeros(1,diff)];

            y = x1_tx + x2_tx;

            if isempty(y_store) && c == 1 && Q == 64 && S == 20
                y_store = y;
            end

            x1_rx = y .* exp(-1i*omega1*(1:length(y)));
            x1_rx = conv(x1_rx, G1);
            x1_rx = downsample(x1_rx, M1) * (1/A1);

            x2_rx = y .* exp(-1i*omega2*(1:length(y)));
            x2_rx = conv(x2_rx, G2);
            x2_rx = downsample(x2_rx, M2) * (1/A2);

            %L = min([length(x1_rx), length(x1)]);
            x1_est = x1_rx(S+1:S+L1);
            

            %L = min([length(x2_rx), length(x2)]);
            x2_est = x2_rx(S+1:S+L2);

            SIDR1(S) = 10*log10(sum(abs(x1).^2) / sum(abs(x1 - x1_est).^2));
            SIDR2(S) = 10*log10(sum(abs(x2).^2) / sum(abs(x2 - x2_est).^2));

            err1(S) = length(find(qamdemod(x1,Q) ~= qamdemod(x1_est,Q)));
            err2(S) = length(find(qamdemod(x2,Q) ~= qamdemod(x2_est,Q)));

        end

        SIDR1_all(c,:) = SIDR1;
        SIDR2_all(c,:) = SIDR2;

        if A1==1 && A2==1, figure(4);
        elseif A1==1 && A2==0.1, figure(6);
        else, figure(8);
        end

        subplot(2,1,1+(Q==64));
        plot(S1,SIDR1); hold on;
        plot(S1,SIDR2);
        title("SIDR A1=" + A1 + " A2=" + A2 + " Q=" + Q);
        legend("x1","x2");

        if A1==1 && A2==1, figure(5);
        elseif A1==1 && A2==0.1, figure(7);
        else, figure(9);
        end

        subplot(2,1,1+(Q==64));
        plot(S1,err1); hold on;
        plot(S1,err2);
        title("Errors A1=" + A1 + " A2=" + A2 + " Q=" + Q);
        legend("x1","x2");

    end
end

NFFT = 2^nextpow2(length(y_store));
Y = fftshift(fft(y_store, NFFT));

fs = 30e6;

f = linspace(-fs/2, fs/2, NFFT);

figure;
plot(f/1e6, 20*log10(abs(Y)));
grid on;

title('Spectrum of y_{tx}(n)');
xlabel('Frequency (MHz)');
ylabel('Magnitude (dB)');



A1 = 1;
S = 20;

A2_vec = logspace(-6,0,120);

figure;
hold on;
set(gca,'XScale','log');

for Q = [4 64]

    x1 = qammod(randi([0 Q-1],L1,1)',Q);
    x2 = qammod(randi([0 Q-1],L2,1)',Q);

    e1 = zeros(size(A2_vec));
    e2 = zeros(size(A2_vec));

    for k = 1:length(A2_vec)

        A2 = A2_vec(k);

        G1 = rcosdesign(rolloff1,S,M1,'sqrt')/sqrt(M1);
        G2 = rcosdesign(rolloff2,S,M2,'sqrt')/sqrt(M2);

        H1 = M1*G1;
        H2 = M2*G2;

        x1_tx = upsample(x1,M1);
        x1_tx = conv(x1_tx,H1)*A1;
        x1_tx = x1_tx .* exp(1i*omega1*(1:length(x1_tx)));

        x2_tx = upsample(x2,M2);
        x2_tx = conv(x2_tx,H2)*A2;
        x2_tx = x2_tx .* exp(1i*omega2*(1:length(x2_tx)));

        diff = length(x2_tx)-length(x1_tx);
        x1_tx = [x1_tx zeros(1,diff)];

        y = x1_tx + x2_tx;

        x1_rx = conv(y.*exp(-1i*omega1*(1:length(y))),G1);
        x1_rx = downsample(x1_rx,M1)/A1;
        x1_est = x1_rx(S+1:end-2*S);

        x2_rx = conv(y.*exp(-1i*omega2*(1:length(y))),G2);
        x2_rx = downsample(x2_rx,M2)/A2;
        x2_est = x2_rx(S+1:end-S);

        L = min([length(x1_est),length(x1)]);
        e1(k) = length(find(qamdemod(x1(1:L),Q) ~= qamdemod(x1_est(1:L),Q)));

        L = min([length(x2_est),length(x2)]);
        e2(k) = length(find(qamdemod(x2(1:L),Q) ~= qamdemod(x2_est(1:L),Q)));

    end

    semilogx(A2_vec,e1);
    semilogx(A2_vec,e2);

end

xlabel("A2");
ylabel("Symbol errors");
title("A1=1, S=20");
legend("x1 Q4","x2 Q4","x1 Q64","x2 Q64");