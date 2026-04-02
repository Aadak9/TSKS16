L=8192;
n=0:L-1;
x=sin(0.25*pi*n)+0.1*sin(0.5*pi*n);
w=window(’blackmanharris’,L)’; %The “’” makes w a row vector like x
xw=x.*w;
wT=linspace(0,2*pi-2*pi/L,L);
Xw=fft(xw);
A=20*log10(abs(Xw)/max(abs(Xw)));
figure;
plot(wT/pi,A);
xlabel(’Normalized frequency, \omegaT/\pi’);
ylabel(’Amplitude spectrum’)