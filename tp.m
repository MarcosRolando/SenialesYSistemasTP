pkg load signal
pkg load matgeom % Para usar la funcion projPointOnLine() y createLine()
pkg load communications % 
global audios Fs c d hist_figure_counter delay_time_counter

c = 340; % Velocidad del sonido en el aire en m/s
d = 0.05; % Distancia entre los microfonos en m
hist_figure_counter = 100; % Para crear distintas figuras para los histogramas
delay_time_counter = 50; % Para crear distintas figuras para los retardos en el tiempo

function plot_lines(slope)
    hold
    grid on

    for i = (1:4)
        line([(2 + 0.05*i) 0], [1 (1 - slope(i)*(2 + 0.05*i))])  
    endfor
    plot([2, 2.05, 2.1, 2.15, 2.2], [1, 1, 1, 1, 1], 'x', 'color', 'r') 
    plot([0.83], [3.40], 'o', 'color', 'r')
    xlim([0, 3])
    ylim([0, 4])
    xlabel("x [m]")
    ylabel("y [m]")
endfunction

%load audios1.mat % Cargo las seniales sin ruido
load audiosRuido1.mat % Cargo las seniales con ruido

close all % Para cerrar los graficos al correr de nuevo el programa

audios(384000:end,:) = []; % Elimino los ultimos 2 segundos aprox (del 8 al 10 en segundos) donde la senial es 0
audios(1:25000,:) = [];

% L = 10;
% b = fir1(150, 1/L, "low") * L; % Filtro pasabajos interpolador para el upsampling
% delay = round(mean(grpdelay(b)));

% for i = (1:5)
%     aux = upsample(audios(:,i), L); % Upsample de gph
%     gph = zeros(rows(aux)+delay, 1); % Agrego delay cantidad de 0s de forma que no pierda la informacion original en el filtro al final de la senial
%     gph(1:end-delay) = aux;
%     gph = filter(b, 1 , gph); % Filtramos/Interpolamos
%     gph = gph(delay+1:end); % Anulamos el desfase introducido por el filtro
%     audios_f(:,i) = gph;
% endfor
% audios = audios_f;

% Ejercicio 5
%{
for i = (1:5)
    audios(:,i) = awgn(audios(:,i), 20, "measured", 1); % Agrego 20dB de ruido blanco a cada senial
endfor
%}

% Ejercicio 1

Fs = 48000; % 48 KHz

t = (0:rows(audios)-1)/Fs;

hf = figure(1);
hold
grid on

for i = (1:5)
    plot(t, audios(:,i) + (i - 1)*0.12)
endfor

xlim([0, 7.45])
ylim([-0.1, 0.75])

legend("Micrófono 1", "Micrófono 2",
        "Micrófono 3", "Micrófono 4"
        , "Micrófono 5")
xlabel("Tiempo [s]")
ylabel("Amplitud")

print(hf, "audios_con_ruido.png")

% Hacemos Zoom para apreciar el retardo

hf = figure(2);
hold
grid on

for i = (1:5)
    %stem(t, audios(:,i))
    plot(t, audios(:, i))
endfor

xlim([4430/Fs, 4451/Fs])
ylim([-0.0024, 0])

legend("Micrófono 1", "Micrófono 2",
        "Micrófono 3", "Micrófono 4"
        , "Micrófono 5")
xlabel("Tiempo [s]")
ylabel("Amplitud")

print(hf, "audios_con_ruido_zoom.png")

% Amplitud
hf = figure(3);
hold
grid on
for i = (1:5)
    fourier_audio = abs(fftshift(fft(audios(:,i))));
    plot((0:1:rows(fourier_audio)-1)*Fs/rows(fourier_audio) - Fs/2, fourier_audio + (i-1)*200)
endfor
xlim([-4000, 4000])
legend("Micrófono 1", "Micrófono 2",
        "Micrófono 3", "Micrófono 4"
        , "Micrófono 5")
xlabel("Frecuencia [Hz]")
ylabel("Amplitud")

print(hf, "transf_fourier_con_ruido.png")

hf = figure(4);
specgram(audios(:,1), 2500, Fs, hanning(2500))
xlabel("Tiempo [s]")
ylabel("Frecuencia [Hz]")
ylim([0, 6500])

print(hf, "specgram_con_ruido.png")

% Ejercicio 2

tau = [-60*10^(-6), -60*10^(-6), -60*10^(-6), -60*10^(-6)]; % Taus medidos visualmente con el grafico del Ejercicio 1
tita = acos(tau * c / d);
slope = tan(tita);

hf = figure(20);
plot_lines(slope);
print(hf, "ejercicio2_sin_ruido.png")

% Ejercicio 3

for i = (1:4)
    [coeffs, lags] = xcorr(audios(:,i+1), audios(:,i), 20); % Paso 
    [~, max_index] = max(coeffs);
    k(i) = -lags(max_index); % Dado que xcorr usa en la definicion x[i + k] pero el enunciado usa x[i - k],
                              % debemos hacer un negado del indice para que sea equivalente a la del enunciado
endfor

format short e

disp("\nRetardos obtenidos mediante el primer método (en segundos):\n")
taus = k / Fs % Retardos entre los microfonos consecutivos

for i = (1:4)
    dft_curr_mic = fft(audios(:,i));
    dft_next_mic = fft(audios(:,i+1));
    R = dft_curr_mic .* conj(dft_next_mic);
    Gph = R ./ abs(R);
    gph = real(ifft(Gph));
    [~, max_index] = max(gph);
    if (max_index > rows(audios)/2)
        max_index = max_index - rows(audios);
    endif
    taus(i) = max_index - 1;
endfor

disp("\nRetardos obtenidos mediante el segundo método (en segundos):\n")
taus = taus / Fs % Paso de muestras a tiempo

% Ejercicio 4

function slope = calculate_lines (audios, Fs, N, delta_n, upsample_gph)
    global c d hist_figure_counter delay_time_counter
    L = 1;

    if (upsample_gph)
        L = 10; % Factor de upsampling
        Fs = Fs * L; % Nueva frecuencia de muestreo
        b = fir1(150, 1/L, "low") * L; % Filtro pasabajos interpolador para el upsampling
        delay = mean(round(grpdelay(b))); % Retardo del filtro, como es lineal el filtro entonces el valor medio es el valor de cada elemento realmente y ademas es la mitad del orden del filtro
        for i = (1:20)
            hist_bins(i) = (i + 29) / Fs;
        endfor
    else
        for i = (1:5)
            hist_bins(i) = i / Fs;
        endfor
    endif

    tita = []; % Angulos
    slope = [];% Pendientes

    for i = (1:4)
        x = audios(:,i);
        y = audios(:,i+1);
        n0 = N / 2;
        j = 1;

        while ((n0 + N/2) < rows(audios))
            w_start = n0 - N/2 + 1;
            w_end = n0 + N/2;
            dft_x = fft(x(w_start:w_end));
            dft_y = fft(y(w_start:w_end));
            R = dft_x .* conj(dft_y);
            Gph = R ./ abs(R); % Un epsilon que evita la inestabilidad del algoritmo cuando R es muy cercano a 0
            gph = real(ifft(Gph));
            if (upsample_gph) % Realizar un upsample de gph para estimar mejor los retardos
                aux = upsample(gph, L); % Upsample de gph
                gph = zeros(rows(aux)+delay, 1); % Agrego delay cantidad de 0s de forma que no pierda la informacion original en el filtro al final de la senial
                gph(1:end-delay) = aux;
                gph = filter(b, 1 , gph); % Filtramos/Interpolamos
                gph = gph(delay+1:end); % Anulamos el desfase introducido por el filtro
            endif
            [~, m] = max(gph);
            if (m > N*L/2) % Si el indice es mayor a la mitad del ancho de la ventana
                m = m - N*L;
            endif
            curr_tau = (m-1) / Fs;
            if (m < 1) % No tiene sentido que me de retardo inverso o que no hay retardo asi que no lo considero
                tau_xy(j) = (m-1) / Fs;
                j = j + 1;
            endif
            n0 = n0 + delta_n; % Actualizo la ventana
        endwhile

        hf = figure(hist_figure_counter);
        h = gca();
        hist(abs(tau_xy), hist_bins) % Grafico el histograma de los retardos
        xlabel(h, "Retardos [s]")
        ylabel(h, "Cantidad")
        print(hf, strcat("hist_audios_con_ruido", num2str(hist_figure_counter), ".png"))
        hf = figure(delay_time_counter);
        hold
        h = gca();
        plot(abs(tau_xy), 'o')
        mean_tau = mean(tau_xy)
        line([0 length(tau_xy)], abs([mean_tau mean_tau]), 'color', 'r')
        ylabel(h, "Retardo [s]")
        print(hf, strcat("delay_audios_con_ruido", num2str(delay_time_counter), ".png"))
        delay_time_counter = delay_time_counter + 1;
        hist_figure_counter = hist_figure_counter + 1;
        tita(i) = acos(mean_tau * c / d);
        slope(i) = tan(tita(i));
        clear tau_xy
    endfor

endfunction

N = 20000;
delta_n = N/2;
slope = calculate_lines(audios, Fs, N, delta_n, false);

hf = figure(5);
plot_lines(slope);
print(hf, "rectas_ejercicio4_con_ruido_con_fuente.png")

%Ejercicio 6

% disp("\n\nRetardos ejercicio 6")
% slope = calculate_lines(audios, Fs, N, delta_n, true);

% hf = figure(6);
% plot_lines(slope)
% print(hf, "rectas_ejercicio4_con_ruido_con_fuente_con_sobremuestreo.png")

% Ejercicio 7

b = fir1(150, 5000/(Fs / 2), "low"); % Disenio el filtro pasabanda que filtre desde los 80Hz hasta los 800Hz
delay = round(mean(grpdelay(b))); % Es la mitad del orden del filtro
audios_noise_filtered = [];
for i = (1:5)
    aux = audios(:,i);
    aux(end+1:end+delay) = 0;
    audios_noise_filtered(:,i) = filter(b, 1, aux);
endfor
audios_noise_filtered = audios_noise_filtered(delay+1:end, :);

% hf = figure(7);
% specgram(audios_noise_filtered(:,1), 2500, Fs, hanning(2500))
% ylim([0, 6500])
% print(hf, "ejercicio7_specgram.png")

% hf = figure(8);
% hold
% freqz(b, 1, length(b), Fs);
% print(hf, "ejercicio7_respuesta_filtro.png")

% hf = figure(9);
% zplane(roots(b));
% print(hf, "ejercicio7_polos_y_ceros.png")

% disp("\nRetardos ejercicio 7")
% slope = calculate_lines(audios_noise_filtered, Fs, N, delta_n, true);
% hf = figure(10);
% plot_lines(slope)
% print(hf, "ejercicio7_posicion_fuente.png")

% Upsampling de las seniales originales

L = 10;
b = fir1(150, 1/L, "low") * L; % Filtro pasabajos interpolador para el upsampling
delay = round(mean(grpdelay(b)));

for i = (1:5)
    aux = upsample(audios_noise_filtered(:,i), L); % Upsample de gph
    gph = zeros(rows(aux)+delay, 1); % Agrego delay cantidad de 0s de forma que no pierda la informacion original en el filtro al final de la senial
    gph(1:end-delay) = aux;
    gph = filter(b, 1 , gph); % Filtramos/Interpolamos
    gph = gph(delay+1:end); % Anulamos el desfase introducido por el filtro
    audios_f(:,i) = gph;
endfor

for i = (1:4)
    [coeffs, lags] = xcorr(audios_f(:,i+1), audios_f(:,i), 100); % Paso 
    [~, max_index] = max(coeffs);
    k(i) = -lags(max_index); % Dado que xcorr usa en la definicion x[i + k] pero el enunciado usa x[i - k],
                              % debemos hacer un negado del indice para que sea equivalente a la del enunciado
endfor

format short e

Fs = Fs * L;

disp("\nRetardos obtenidos mediante el primer método (en segundos):\n")
taus = k / Fs % Retardos entre los microfonos consecutivos

hf = figure(11);
tita = acos(taus * c / d);
slope = tan(tita);
plot_lines(slope)
print(hf, "ejercicio7_mejor_fuente.png")

% Ejercicio 8

aux = zeros(rows(audios_f), 5);
aux(:,1) = audios_f(:,1);
delay = 0;
for i = (2:5)
    delay = delay + k(i-1);
    aux(1:end-abs(delay),i) = audios_f(1+abs(delay):end,i);
endfor
audios_f = aux;
audio_prom = audios_f(:,1);
for i = (2:5)
    audio_prom = audio_prom + audios_f(:,i); 
endfor
audio_prom = audio_prom / 5;


clear all % Clear all variables
