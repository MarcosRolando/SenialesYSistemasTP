pkg load signal
pkg load matgeom % Para usar la funcion projPointOnLine() y createLine()
pkg load communications % 
global audios Fs

%load audios1.mat % Cargo las seniales sin ruido
load audiosRuido1.mat % Cargo las seniales con ruido

close all % Para cerrar los graficos al correr de nuevo el programa

audios(400000:end,:) = []; % Elimino los ultimos 2 segundos aprox (del 8 al 10 en segundos) donde no hay senial 

% Ejercicio 5
%{
for i = (1:5)
    audios(:,i) = awgn(audios(:,i), 20, "measured", 1); % Agrego 20dB de ruido blanco a cada senial
endfor
%}

% Ejercicio 1

Fs = 48 * 10^3; % 48 KHz

t = (0:rows(audios)-1)/Fs;

figure(1)
hold
grid on

for i = (1:5)
    plot(t, audios(:,i) + (i - 1)*0.1)
endfor

xlim([0, 8])

legend("Micrófono 1", "Micrófono 2",
        "Micrófono 3", "Micrófono 4"
        , "Micrófono 5")

% Zoom in to appreciate the offset

figure(2)
hold
grid on

for i = (1:5)
    plot(t, audios(:,i))
endfor

xlim([29430/Fs, 29451/Fs])
ylim([-0.0024, 0])

legend("Micrófono 1", "Micrófono 2",
        "Micrófono 3", "Micrófono 4"
        , "Micrófono 5")


figure(3)
grid on
fourier_audio = abs(fftshift(fft(audios(:,1))));
plot((0:1:rows(fourier_audio)-1)*2*pi/rows(fourier_audio) - pi, fourier_audio)
xlim([-pi/5, pi/5])
legend("Fourier Transform")

figure(4)
specgram(audios(:,1), 2500, Fs, hanning(2500))
xlim([0.4, 8])
ylim([0, 1500])

% Ejercicio 2

% Estimar los retardos viendo el grafico del segmento ese del ejercicio 1

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
    Gph = dft_curr_mic .* conj(dft_next_mic) ./ (abs(dft_curr_mic).*abs(dft_next_mic));
    gph = real(ifft(Gph));
    [~, max_index] = max(gph);
    if (max_index > rows(audios)/2)
        max_index = max_index - rows(audios);
    endif
    taus(i) = max_index;
endfor

disp("\nRetardos obtenidos mediante el segundo método (en segundos):\n")
taus = taus / Fs % Paso de muestras a tiempo

% Ejercicio 4

function slope = calculate_lines (N, delta_n)
    global audios Fs

    c = 340; % Velocidad del sonido en el aire en m/s
    d = 0.05; % Distancia entre los microfonos en m

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
            Gph = dft_x .* conj(dft_y) ./ (abs(dft_x) .* abs(dft_y));
            gph = real(ifft(Gph));
            [~, m] = max(gph);
            if (m > N/2)
                m = m - N;
            endif
            tau_xy(j) = m / Fs;
            j = j + 1;
            n0 = n0 + delta_n; 
        endwhile

        mean_tau = mean(tau_xy)
        tita(i) = acos(mean_tau * c / d);
        slope(i) = tan(tita(i));
        clear tau_xy
    endfor

endfunction


% best_error = 10^4; % Un error grande para que la primera iteracion siempre lo reemplace
% for N = (1000:1000:100000)
%     for delta_n = ((N/10):(N*5/100):N)
%         slope = calculate_lines(N, delta_n);

%         k = 1;
%         best_x = 0;
%         best_y = 10^3; % Un valor grande para que la primera iteracion siempre lo reemplace
%         for i = (1:3)
%             for j = (i:3)
%                 if (slope(i) == slope(j+1))
%                     continue
%                 endif
%                 curr_x = (-slope(j+1)*0.05*(j+1) + slope(i)*0.05*i) / (slope(i) - slope(j+1)); 
%                 curr_y = slope(i)*curr_x - slope(i)*0.05*i;
%                 if (curr_y > 0 && curr_y < best_y)
%                     best_x = curr_x;
%                     best_y = curr_y;
%                 endif
%                 k = k + 1;
%             endfor
%         endfor
        
%         curr_error = 0;
%         for i = (1:4)
%             Line = createLine([0.05*i 0], [-1 (slope(i)*(-1) - slope(i)*0.05*i)]);
%             v = [best_x best_y] - projPointOnLine([best_x best_y], Line);
%             curr_error = curr_error + sqrt(v(1)^2 + v(2)^2);
%         endfor
%         if (curr_error < best_error)
%             best_error = curr_error
%             best_N = N
%             best_delta_n = delta_n
%         endif

%     endfor
% endfor


best_N = 45000; % Calculados ejecutando el algoritmo de arriba
best_delta_n = 33750;
slope = calculate_lines(best_N, best_delta_n);

figure(5)
hold
grid on

for i = (1:4)
    line([0.05*i -1], [0 slope(i) * (-1 - 0.05 * i)])  
endfor

% Ejercicio 6



clear all % Clear all variables
