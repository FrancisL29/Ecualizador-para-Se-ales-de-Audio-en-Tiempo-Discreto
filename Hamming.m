clc; clear; close all;

N = 301;                 % Longitud de la ventana
w = hamming(N);          % Ventana Hamming

n = 0:N-1;

figure
plot(n, w, 'LineWidth',2)
grid on
xlabel('Muestras (n)')
ylabel('Amplitud')
title('Ventana Hamming en el dominio del tiempo')
Nfft = 8192;                 % FFT grande = mejor resolución
W = fft(w, Nfft);            
W = fftshift(W);             % Centrar 0 Hz

f = linspace(-0.5,0.5,Nfft); % Frecuencia normalizada

figure
plot(f, 20*log10(abs(W)/max(abs(W))), 'LineWidth',2)
grid on
xlabel('Frecuencia normalizada (×π rad/muestra)')
ylabel('Magnitud (dB)')
title('Espectro de la ventana Hamming')
ylim([-120 5])

