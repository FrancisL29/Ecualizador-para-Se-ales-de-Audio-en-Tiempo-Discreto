%% =================== ECUALIZADOR 6 BANDAS + PASABAJO SRC
clear; close all; clc;

Fs = 44100;

%% ===== BANDAS DEL ECUALIZADOR
bandsHz = [16 60; 60 250; 250 2000; 2000 4000; 4000 6000; 6000 16000];
numBands = size(bandsHz,1);

%% ===== GANANCIAS POR BANDA (dB)
GdB = [0, 0, 10, 0, 25, 15];   % Modifica entre -25 y 25
Glin = 10.^(GdB/20);

%% ===== DISEÑO DE FILTROS PASABANDA
hBP = cell(numBands,1);
for k = 1:numBands
    hBP{k} = bandpass_fir(bandsHz(k,1), bandsHz(k,2), Fs);
end

%% ===== FILTRO PASABAJO (SRC / ANTIALIAS)
factor = 4;
Lf = 101;
n = -(Lf-1)/2:(Lf-1)/2;
wc = pi/factor;
hLP = (wc/pi)*sinc((wc/pi)*n).*hamming(Lf)';

%% ===== SEÑAL MULTIBANDA DE PRUEBA
N = 16384;
t = (0:N-1)/Fs;
x = sin(2*pi*40*t) + sin(2*pi*120*t) + sin(2*pi*2000*t) + ...
    sin(2*pi*4000*t) + sin(2*pi*6000*t) + sin(2*pi*10000*t);

f = (0:N-1)*(Fs/N);

%% ===== SEÑAL ORIGINAL EN FRECUENCIA
X = abs(fft(x))/N;
figure
plot(f(1:N/2),20*log10(X(1:N/2)),'LineWidth',1.5)
grid on
xlabel('Frecuencia (Hz)')
ylabel('Magnitud (dB)')
title('Señal Original – Dominio de Frecuencia')
xlim([0 16000])

%% ===== ESPECTRO DE FILTRO PASABANDA (CON GANANCIA)
figure
hold on
for k = 1:numBands
    [H,freq] = freqz(hBP{k},1,4096,Fs);
    plot(freq,20*log10(abs(H)*Glin(k)),'LineWidth',1.5)
end
grid on
xlabel('Frecuencia (Hz)')
ylabel('Magnitud (dB)')
title('Espectro de Filtros Pasabanda con Ganancia')
legend('16–60','60–250','250–2000','2000–4000','4000–6000','6000–16000','Location','Best')

%% ===== SEÑAL DESPUÉS DE PASAR POR EL ECUALIZADOR
yEQ = zeros(size(x));
for k = 1:numBands
    yEQ = yEQ + filter(Glin(k)*hBP{k},1,x);
end
Yeq = abs(fft(yEQ))/N;

figure
plot(f(1:N/2),20*log10(Yeq(1:N/2)),'LineWidth',1.5)
grid on
xlabel('Frecuencia (Hz)')
ylabel('Magnitud (dB)')
title('Señal después de pasar por Filtros Pasabanda')
xlim([0 16000])

%% ===== ️FILTRO PASABAJO
figure
[Hlp,freq] = freqz(hLP,1,4096,Fs);
plot(freq,20*log10(abs(Hlp)),'k','LineWidth',2)
grid on
xlabel('Frecuencia (Hz)')
ylabel('Magnitud (dB)')
title('Filtro Pasabajo')
xlim([0 16000])

%% ===== SEÑAL FINAL DESPUÉS DEL PASABAJO
yLP = filter(hLP,1,yEQ);
Ylp = abs(fft(yLP))/N;

figure
plot(f(1:N/2),20*log10(Ylp(1:N/2)),'LineWidth',1.5)
grid on
xlabel('Frecuencia (Hz)')
ylabel('Magnitud (dB)')
title('Señal Final después del Filtro Pasabajo')
xlim([0 16000])


%% =RESPUESTA AL IMPULSO

figure;
for k = 1:numBands
    subplot(numBands,1,k);
    stem(hBP{k},'filled');
    title(['Respuesta al impulso del filtro pasabanda ' num2str(k)]);
    xlabel('Muestras');
    ylabel('Amplitud');
end

figure;
stem(hLP,'filled');
title('Respuesta al impulso del filtro pasabajo');
xlabel('Muestras');
ylabel('Amplitud');



%% ===== RESPUESTA TOTAL DEL ECUALIZADOR (IMPULSO)
hEQ = zeros(size(hBP{1}));

for k = 1:numBands
    hEQ = hEQ + Glin(k)*hBP{k};
end

figure;
stem(hEQ,'filled')
title('Respuesta al impulso TOTAL del Ecualizador FIR')
xlabel('Muestras')
ylabel('Amplitud')
grid on


figure;
[Heq,freq] = freqz(hEQ,1,4096,Fs);
plot(freq,20*log10(abs(Heq)),'LineWidth',2)
grid on
xlabel('Frecuencia (Hz)')
ylabel('Magnitud (dB)')
title('Respuesta en Frecuencia TOTAL del Ecualizador')
xlim([0 16000])





%% ===== FUNCION FILTRO PASABANDA
function h = bandpass_fir(f1,f2,Fs)
    L = 201;
    n = -(L-1)/2:(L-1)/2;
    w1 = 2*pi*f1/Fs;
    w2 = 2*pi*f2/Fs;
    h = (w2/pi)*sinc((w2/pi)*n) - (w1/pi)*sinc((w1/pi)*n);
    h = h .* hamming(L)';
end


