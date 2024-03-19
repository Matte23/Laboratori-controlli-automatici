%% Setup dati
s = tf('s');

A = [0 -1 5; 0 0 3; 0 0 -2];
B = [1 1 1]';
C = [0 0 5];
D = 0;

W = minreal(zpk(inv(s*eye(size(A))-A)));
% Nota: valutare sempre che faccia le cancellazioni
% se no inserire una tolleranza

% Definizione t
t = linspace(0,4,100);

%% Punto 1
%%% Stabilità interna
E = real(eig(A))
% Il -2 appare con molteplicità singola -> OK
% Lo zero appare con molteplicità doppia, dobbiamo valutare il polinomio minimo

% Calcolo diretto delle radici del polinomio minimo con MATLAB
E_min = roots(minpoly(A))

% Alternativa calcolando il polinomio minimo identificando il lcm della matrice Af
% Per semplificare numeratore / denominatore prima fattorizzo con zpk e poi chiamo minreal
Af = W;
% il polinomio minimo è
p_min = s^2*(s+2);

% Lo 0 appare con molteplicità due, pertanto il sistema è instabile

%%% Stabilità BIBO
H = C*W*B;
p = pole(H)

% Tutti i poli devono essere negativi
% L'unico polo è -2 (negativo), quindi è BIBO stabile

%% Punto 2
% Calcolo risposta forzata con u(t) = 9 epsilon(t)
U = 9/s; % Trasformata dell'ingresso
Yf = minreal(zpk(H*U));

[num_Xf,den_Xf] = tfdata(Yf,'v');
[r1,p1] = residue(num_Xf,den_Xf)
% Ys = 22.5/s - 22.5/(s+2)
Ytf = 22.5 - 22.5 * exp(-2*t);

%% Punto 3
% Calcolo risposta libera (dati 1)
Xz1 = [1 5 0]';
Yf = minreal(zpk(C*W*Xz1));

[num_Xf,den_Xf] = tfdata(Yf,'v');
[r1,p1] = residue(num_Xf,den_Xf)
Yt1 = 0*t;
% Perché X(3) dipende solo da se stesso e dall'ingresso (e non da X(1) e X(2))

%% Punto 4
% Calcolo risposta libera (dati 2)
Xz1 = [0 0 3]';
Yf = minreal(zpk(C*W*Xz1));

[num_Xf,den_Xf] = tfdata(Yf,'v');
[r1,p1] = residue(num_Xf,den_Xf)
Yt2 = 15 * exp(-2*t);

%% Punto 5
% Calcoli dei punti 2 3 4 usando ss()

% Calcolo risposta forzata con u(t) = 9 epsilon(t)
S = ss(A,B,C,D);
% Nota: possiamo anche calcolare la funzione di trasferimento con
H1 = tf(S);
H2 = zpk(S);
Yf = step(9*S,t);

% Calcolo risposta libera (dati 1)
Xi1 = [1 5 0]';
Y1 = initial(S, Xi1, t);

% Calcolo risposta libera (dati 2)
Xi2 = [0 0 3]';
Y2 = initial(S, Xi2, t);

% Plot risposta forzata
figure
hold on;
plot(t,Yf, 'b');
plot(t,Ytf, '.r');
title("Risposta forzata")
xlabel('Tempo (secondi)')
ylabel('Ampiezza')
grid

% Plot due risposte libere
figure
tiledlayout(2,1)

nexttile
hold on
plot(t,Y1, 'b');
plot(t,Yt1, '.r');
title("Risposta libera X(0) = [1 5 0]'")
xlabel('Tempo (secondi)')
ylabel('Ampiezza')
grid

nexttile;
hold on
plot(t,Y2, 'b');
plot(t,Yt2, '.r');
title("Risposta libera X(0) = [0 0 3]'")
xlabel('Tempo (secondi)')
ylabel('Ampiezza')
grid

%% Punto bonus: somma risposta forzata con risposte libere
figure
tiledlayout(2,1)

nexttile
plot(t,Y1+Yf);
title("Risposta libera X(0) = [1 5 0]' con risposta forzata")
xlabel('Tempo (secondi)')
ylabel('Ampiezza')
grid

nexttile;
plot(t,Y2+Yf);
title("Risposta libera X(0) = [0 0 3]' con risposta forzata")
xlabel('Tempo (secondi)')
ylabel('Ampiezza')
grid
