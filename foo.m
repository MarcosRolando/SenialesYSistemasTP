load audios1.mat % Loads the file with the signals

figure(1)
hold
grid on

for i = (1:5)
    plot(audios(:,i) + (i - 1)*0.1)
endfor

% Zoom in to appreciate the offset

figure(2)
hold
grid on

for i = (1:5)
    plot(audios(:,i))
endfor

xlim([29430, 29451])
ylim([-0.0024, 0])

