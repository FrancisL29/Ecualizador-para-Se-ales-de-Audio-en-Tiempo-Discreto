function terminado_comen           
clc; close all;                    

%% ======================== FIGURA =======================
fig = figure('Name','DSP Audio System','Position',[30 30 1200 700]);
% Crea la ventana principal del programa (interfaz gráfica)

axTime = axes(fig,'Position',[0.05 0.55 0.35 0.35]);
% Crea un eje donde se dibujará la señal en el tiempo
title(axTime,'x[n] - Dominio del tiempo'); grid on;

axSpecOrig = axes(fig,'Position',[0.05 0.08 0.35 0.30]);
% Eje donde se mostrará el espectro de la señal original
title(axSpecOrig,'|X(e^{j\omega})| original'); grid on;

axSpecProc = axes(fig,'Position',[0.45 0.08 0.35 0.30]);
% Eje donde se mostrará el espectro de la señal procesada
title(axSpecProc,'|Y(e^{j\omega})| procesado'); grid on;




%% ============================ CONTROLES ==============================
uicontrol(fig,'Style','pushbutton','String','Cargar Audio',...
    'Position',[500 600 130 40],'Callback',@CargarAudio);


modeMenu = uicontrol(fig,'Style','popupmenu',...
    'String',{'Decimación','Interpolación'},...
    'Position',[500 540 150 30],'Callback',@CambiarModo);


uicontrol(fig,'Style','text','String','Factor M / L',...
    'Position',[500 510 120 20]);

sliderML = uicontrol(fig,'Style','slider','Min',1,'Max',40,'Value',2,...
    'Position',[500 490 150 20],'Callback',@actTODO);

txtML = uicontrol(fig,'Style','text','String','2',...
    'Position',[660 490 40 20]);

uicontrol(fig,'Style','pushbutton','BackgroundColor',[0.2 0.8 0.2],'String','PLAY',...
    'Position',[500 440 120 40],'Callback',@InicioBloques);

uicontrol(fig,'Style','pushbutton', 'BackgroundColor',[1 0 0],'String','STOP',...
    'Position',[630 440 120 40],'Callback',@stopDSP);

uicontrol(fig,'Style','pushbutton','String','RESET',...
    'Position',[500 390 250 40],'Callback',@resetEQ);







%% ============================= ECUALIZADOR ===============================
TextoBandas = {'Sub (16–60)','Bass (60–250)','LowMid (250–2k)',...
            'HighMid (2–4k)','Presence (4–6k)','Brill (6–16k)'};
% Nombres de las 6 bandas de frecuencia

GananSliders = gobjects(6,1);  % Vector para guardar sliders de ganancia
GananTxt     = gobjects(6,1);  % Vector para mostrar valor de ganancia

for k=1:6
    y = 660-40*k;  % Posición vertical de cada control

    uicontrol(fig,'Style','text','String',TextoBandas{k},...
        'Position',[800 y 130 20],'HorizontalAlignment','left');
    % Texto del nombre de la banda

    GananSliders(k) = uicontrol(fig,'Style','slider',...
        'Min',-20,'Max',20,'Value',0,...
        'Position',[940 y 120 20],'Callback',@actTODO);
    % Slider de ganancia en dB

    GananTxt(k) = uicontrol(fig,'Style','text','String','0 dB',...
        'Position',[1080 y 60 20]);
    % Texto que muestra el valor actual de la ganancia
end








%% ================================ VARIABLES ==========================
x = []; Fs = 0;               
EQfiltro = {};               % Aquí se guardarán los filtros del ecualizador
bandaHz = [16 60;60 250;250 2000;2000 4000;4000 6000;6000 16000];


SRFfiltro = [];               % Filtro anti-alias o de interpolación
lastFactor = -1; lastMode = -1; % Guardan valores anteriores para no recalcular filtros
isPlaying = false;            % Indica si el audio se está reproduciendo








%% ================================ CARGA AUDIO ==================================
function CargarAudio(~,~)

    % x = señal (vector), Fs = frecuencia de muestreo
    [x,Fs] = audioread('locos.wav');

    x = mean(x,2);

    % Limita la duración 
    x = x(1:min(end,30*Fs));

    % Vector de tiempo para graficar (en segundos)
    t = (0:length(x)-1)/Fs;

    % Grafica la señal en el dominio del tiempo
    plot(axTime,t,x);

    % FFT centrada (magnitud)
    X = fftshift(abs(fft(x)));

    % Eje de frecuencia centrado desde -Fs/2 a Fs/2
    f = linspace(-Fs/2,Fs/2,length(X));

    % Grafica el espectro original
    plot(axSpecOrig,f,X); xlabel('t');

    % Precalcula los filtros del ecualizador (6 bandas)
    EQfiltro = cell(6,1);
    for b=1:6
        EQfiltro{b} = bandpass_fir(bandaHz(b,1),bandaHz(b,2),Fs);
    end

    % Actualiza todo (ecualización, decimación/interpolación, gráficas)
    actTODO

    % --- Mostrar en consola ---
    fprintf('Frecuencia de muestreo (Fs): %.2f Hz\n', Fs);
    fprintf('Frecuencia máxima (Fs/2): %.2f Hz\n', Fs/2);

end









%% ============================ FILTRO PASABANDA ==========================================
function h = bandpass_fir(f1,f2,Fs)

    L = 301;                         % Longitud del filtro FIR (número de coeficientes).
    n = -(L-1)/2:(L-1)/2;            % Vector de índices centrado (simetría alrededor de cero).
    
    w1 = 2*pi*f1/Fs; 
    w2 = 2*pi*f2/Fs; 

    h = (w2/pi)*sinc((w2/pi)*n) - (w1/pi)*sinc((w1/pi)*n);
    % Construye la respuesta ideal del filtro pasabanda:
    % - Primer término: filtro paso bajo con corte en f2
    % - Segundo término: filtro paso bajo con corte en f1
    % La resta genera un filtro pasabanda (f2 - f1).

    h = h .* hamming(L)';             
    % Aplica ventana Hamming para suavizar el filtro y reducir los
    % lóbulos laterales (reducción de ripple).
end







%% ============================ FILTRO SRC (ANTI-ALIAS / ANTI-IMAGING) =========================
function h = SRC_fir(factor, Lf)

    n = -(Lf-1)/2:(Lf-1)/2;      % Vector de índices centrado
    wc = pi/factor;              % Frecuencia de corte normalizada
    h = (wc/pi)*sinc((wc/pi)*n).*hamming(Lf)';

end








%% ===================================== APLICAR EQ ===================================
function xeq = applyEQ(x)

    xeq = zeros(size(x));           
    % Inicializa la señal resultante con ceros (misma longitud que x).

    for b = 1:6
        g = get(GananSliders(b),'Value'); 
        % Lee la ganancia del slider de la banda b (en dB).

        set(GananTxt(b),'String',[num2str(round(g)) ' dB'])   % Actualiza el texto para mostrar la ganancia actual.
      
        G = 10^(g/20);        % Convierte ganancia de dB a ganancia lineal:       
        
          xeq = xeq + G*conv(x,EQfiltro{b},'same');
        % Aplica el filtro de la banda b mediante convolución y suma el resultado
        % a la señal final (suma de todas las bandas con su ganancia).
    end
end






%% ================================= GRAFICAR TODO =================================
function actTODO(~,~)

    if isempty(x), return; end      % Si no hay señal cargada (x vacío), salir de la función

    % Actualiza el texto del GUI con el valor actual del slider (M o L)
    set(txtML,'String',num2str(round(get(sliderML,'Value'))));

    % Aplica la ecualización a la señal original (x) y guarda el resultado en xeq
    xeq = applyEQ(x);

    % Obtiene el modo seleccionado en el menú (1=Decimación, 2=Interpolación)
    modeVal = get(modeMenu,'Value');

    % Lee el valor del slider (M o L) y lo redondea a entero
    factor  = round(get(sliderML,'Value'));

    % Si el factor o el modo cambiaron desde la última vez:
    if factor~=lastFactor || modeVal~=lastMode
        SRFfiltro = SRC_fir(factor,101);
        % Guarda el factor y modo actuales para evitar recalcular el filtro
        lastFactor = factor; 
        lastMode = modeVal;
    end

    % ----------------- DECIMACIÓN -----------------
    if modeVal==1
        y = conv(xeq,SRFfiltro,'same');  % Filtra la señal para evitar aliasing
        y = y(1:factor:end);             % Reduce la señal tomando 1 muestra cada 'factor'
                                         
        Fs2 = Fs/factor;                 % Nueva frecuencia de muestreo

    % ----------------- INTERPOLACIÓN -----------------
    else
        xu = zeros(length(xeq)*factor,1); % Crea un vector con espacio para la interpolación
        xu(1:factor:end) = xeq;           % Inserta ceros entre muestras (↑L, upsampling)
        y = conv(xu,SRFfiltro,'same');    % Filtra para reconstruir la señal interpolada
        Fs2 = Fs*factor;                  % Nueva frecuencia de muestreo
    end

    
    Y = fftshift(abs(fft(y)));

    % Genera el eje de frecuencia usando la nueva Fs2
    f = linspace(-Fs2/2,Fs2/2,length(Y));

    % Grafica el espectro procesado en el eje correspondiente
    plot(axSpecProc,f,Y); xlabel('Hz');

end





%% ============================== BLOQUES ====================================
function InicioBloques(~,~)

    if isempty(x), return; end
    % Si no señal, salir de la función

    blockSize = 4*Fs;               

    isPlaying = true; % Variable para detener el procesamiento (STOP).
    

    % Recorre la señal en bloques completos (sin considerar el último bloque incompleto)
    for k = 1:floor(length(x)/blockSize)

        fprintf('Numero de bloque: %d\n', k);
        

        if ~isPlaying, break; end % STOP.
        

        xb = x((k-1)*blockSize+1:k*blockSize); % Extrae el bloque k-ésimo de la señal original.
        xb = applyEQ(xb);  % Aplica la ecualización al bloque.
       

        % Lee el modo actual (1=Decimación, 2=Interpolación).
        factor = round(get(sliderML,'Value'));
        modeVal = get(modeMenu,'Value'); 
        

        
        h = SRC_fir(factor,401);

        

        % ---------- DECIMACIÓN ----------
        if modeVal == 1
            yb = conv(xb,h,'same');             % Filtra el bloque para evitar aliasing
            yb = yb(1:factor:end);             % Reduce la tasa de muestreo tomando 1 muestra cada 'factor'.
            Fs2 = Fs/factor;


        % ---------- INTERPOLACIÓN ----------
        else
            xu = zeros(length(xb)*factor,1);            % Crea un vector con espacio para insertar ceros.
            xu(1:factor:end) = xb;                      % Inserta ceros entre muestras (upsampling).
            yb = conv(xu,h,'same');                     % Filtra el resultado para reconstruir la señal interpolada.
            Fs2 = Fs*factor;                            % Nueva frecuencia de muestreo
             
        end

        sound(yb,Fs2); % Reproduce el bloque procesado 
        
        

        %%%%%%%%%%%%%%%% GRAFICAR EN TIEMPO%%%%%%%%%%%%%%%%%%%%%%%%%
%         t = (0:length(yb)-1)/Fs2;   % Eje de tiempo en segundos
%         
%         figure;
%         plot(t, yb);
%         grid on;
%         xlabel('Tiempo [s]');
%         ylabel('Amplitud');
%         title('Señal y_b después de la conversión de tasa de muestreo');
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

          pause(blockSize/Fs);  % Espera el tiempo real equivalente al bloque original
        
 
        
    end
end






%% ================================= STOP ===========================
function stopDSP(~,~)
    isPlaying = false;              % Detiene el bucle
    clear sound                     % Para el audio inmediatamente
end






%% ================================ RESET EQ =============================
function resetEQ(~,~)
    for b = 1:6
        set(GananSliders(b),'Value',0);
        set(GananTxt(b),'String','0 dB');
    end
    actTODO;                      % Recalcula sistema
end




%% =============================== LIMITAR M / L ============================
function CambiarModo(~,~)
    modeVal = get(modeMenu,'Value');

    if modeVal == 1                 % Decimación
        set(sliderML,'Min',1,'Max',40);
        if get(sliderML,'Value') > 40
            set(sliderML,'Value',40);
        end
    else                            % Interpolación
        set(sliderML,'Min',1,'Max',4);
        if get(sliderML,'Value') > 4
            set(sliderML,'Value',4);
        end
    end

    set(txtML,'String',num2str(round(get(sliderML,'Value'))));
    actTODO;
end




end

