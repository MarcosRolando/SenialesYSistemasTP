pkg load signal
load audios1.mat % Loads the file with the signals

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
    [coeffs, lags] = xcorr(audios(:,i), audios(:,i+1), 20);
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
    [~, max_index] = max(gph(1:100)); % Lo limito a 100 puntos porque el maximo esta dentro de ese rango y si tomo todo
                                      % agarraria el ruido del final
    taus(i) = max_index;
endfor

disp("\nRetardos obtenidos mediante el segundo método (en segundos):\n")
taus = taus / Fs % Paso de muestras a tiempo