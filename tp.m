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
    plot([1.56], [2.04], 'o', 'color', 'r')
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

figure(1)
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
xlabel("Tiempo (s)")
ylabel("Amplitud")

% Hacemos Zoom para apreciar el retardo

figure(2)
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
xlabel("Tiempo (s)")
ylabel("Amplitud")

% Amplitud
figure(3)
hold
grid on
for i = (1:5)
    fourier_audio = abs(fftshift(fft(audios(:,i))));
    plot((0:1:rows(fourier_audio)-1)*2*pi/rows(fourier_audio) - pi, fourier_audio + (i-1)*200)
endfor
xlim([-pi/8, pi/8])
legend("Micrófono 1", "Micrófono 2",
        "Micrófono 3", "Micrófono 4"
        , "Micrófono 5")
xlabel("Frecuencia (rad)")
ylabel("Amplitud")

figure(4)
specgram(audios(:,1), 2500, Fs, hanning(2500))
ylim([0, 3000])

% Ejercicio 2

tau = [-60*10^(-6), -60*10^(-6), -60*10^(-6), -60*10^(-6)]; % Taus medidos visualmente con el grafico del Ejercicio 1
tita = acos(tau * c / d);
slope = tan(tita);

figure(20)
plot_lines(slope);

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
            hist_bins(i) = -(i + 29) / Fs;
        endfor
    else
        for i = (1:5)
            hist_bins(i) = -i / Fs;
        endfor
    endif

    hist_bins = flip(hist_bins); % Ordeno el input porque sino me queda al reves
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
            Gph = R ./ abs(R);
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
            if (m < 1) % No tiene sentido que me de retardo inverso o que no hay retardo asi que no lo considero
                tau_xy(j) = (m-1) / Fs;
                j = j + 1;
            endif
            n0 = n0 + delta_n; % Actualizo la ventana
        endwhile

        % figure(hist_figure_counter);
        % h = gca();
        % hist(tau_xy, hist_bins) % Grafico el histograma de los retardos
        % xlabel(h, "Retardos [s]")
        % ylabel(h, "Cantidad")
        % figure(delay_time_counter)
        % hold
        % h = gca();
        % plot(tau_xy, 'o')
        mean_tau = mean(tau_xy)
        % line([0 length(tau_xy)], [mean_tau mean_tau], 'color', 'r')
        % ylabel(h, "Retardo [s]")
        % delay_time_counter = delay_time_counter + 1;
        % hist_figure_counter = hist_figure_counter + 1;
        tita(i) = acos(mean_tau * c / d);
        slope(i) = tan(tita(i));
        clear tau_xy
    endfor

endfunction

N = 20000;
delta_n = N/2;
slope = calculate_lines(audios, Fs, N, delta_n, false);

figure(5)
plot_lines(slope);

% Ejercicio 6

disp("\n\nRetardos ejercicio 6")
slope = calculate_lines(audios, Fs, N, delta_n, true);

figure(6)
plot_lines(slope)

% Ejercicio 7

b = fir1(150, 4000/(Fs / 2), "low"); % Disenio el filtro pasabanda que filtre desde los 80Hz hasta los 800Hz
delay = round(mean(grpdelay(b))); % Es la mitad del orden del filtro
audios_noise_filtered = [];
for i = (1:5)
    aux = audios(:,i);
    aux(end+1:end+delay) = 0;
    audios_noise_filtered(:,i) = filter(b, 1, aux);
endfor
audios_noise_filtered = audios_noise_filtered(delay+1:end, :);

figure(15)
hold
%stem(audios_noise_filtered(end-40:end,1))
%stem(audios(end-40:end,1))
specgram(audios_noise_filtered(:,1), 2500, Fs, hanning(2500))
ylim([0, 10000])

disp("\nRetardos")
slope = calculate_lines(audios_noise_filtered, Fs, N, delta_n, true);

% figure(7)
% hold
% grid on

% for i = (1:4)
%     line([0.05*i -1], [0 slope(i) * (-1 - 0.05 * i)])  
% endfor

clear all % Clear all variables
