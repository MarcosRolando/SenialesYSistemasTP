pkg load signal

close all

global c d

c = 340; % Velocidad del sonido en el aire en m/s
d = 0.20; % Distancia entre los oidos en m

format short e

function plot_lines(slope)
    hold
    grid on
    
    if (slope < 0) % Asi las recta la dibujo en el intervalo de interes
        line([(2 + 0.20) 0], [1 (1 - slope*(2 + 0.20))])  
    else
        line([(2 + 0.20) 4], [1 (1 + slope*(4 - 0.20))])  
    endif
    plot([2, 2.20], [1, 1], 'x', 'color', 'r') 
    xlim([0, 3])
    ylim([0, 4])
    xlabel("x [m]")
    ylabel("y [m]")
endfunction

function slope = calculate_lines (audios, Fs, N, delta_n, upsample_gph)
    global c d
    L = 1;

    if (upsample_gph)
        L = 10; % Factor de upsampling
        Fs = Fs * L; % Nueva frecuencia de muestreo
        b = fir1(150, 1/L, "low") * L; % Filtro pasabajos interpolador para el upsampling
        delay = mean(round(grpdelay(b))); % Retardo del filtro, como es lineal el filtro entonces el valor medio es el valor de cada elemento realmente y ademas es la mitad del orden del filtro
    endif

    tita = []; % Angulos
    slope = [];% Pendientes

    x = audios(:,1);
    y = audios(:,2);
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
        if (m != 1)
            curr_tau = (m-1) / Fs;
            tau_xy(j) = (m-1) / Fs;
            j = j + 1;
        endif
        n0 = n0 + delta_n; % Actualizo la ventana
    endwhile

    mean_tau = mean(tau_xy)
    tita = acos(mean_tau * c / d);
    slope = tan(tita);
    clear tau_xy

endfunction

[audios, Fs] = audioread("audio2.wav");
cut_n = rows(audios) / 5;
N = 20000;
delta_n = N / 2;
disp("\nRetardos medios obtenidos mediante GCC-PHAT (en segundos)\n")
for i = (1:5)
    hf = figure(i);
    slope = calculate_lines(audios(1+cut_n*(i-1):cut_n*i, :), Fs, N, delta_n, true);
    plot_lines(slope)
    print(hf, strcat("ejercicio11", num2str(i), ".png"))
endfor

clear all