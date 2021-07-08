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