%% Setup dati
s = tf('s');

A = [0 -1 5; 0 0 3; 0 0 -2];
B = [1 1 1]';
C = [0 0 5];

Al = inv(s*eye(3)-A);

%% Punto 1
%%% Stabilità interna
E = real(eig(A))
% Lo zero appare con molteplicità doppia, dobbiamo valutare il polinomio minimo

% Calcolo diretto con matlab
roots(minpoly(A))

% Alternativa calcolando a mano il mcm (ovvero il polinomio caratteristico)
% Al = minreal(zpk(Al));

% Lo 0 appare con molteplicità due, pertanto il sistema è instabile

%%% Stabilità BIBO
H = C*Al*B;
p = pole(H)

% Il polo è -2 (negativo), quindi è stabile BIBO

%% Punto 2
% Calcolo risposta forzata con u(t) = 9 g(t)
U = 9/s;
Yf = minreal(zpk(C*Al*B*U));
% In alternativa Yf = H*U;

[num_Xf,den_Xf]=tfdata(Yf,'v');
[r1,p1]=residue(num_Xf,den_Xf)
% Yt = 22.5 - 22.5 * e^(-2t)

%% Punto 3
% Calcolo risposta libera (dati 1)
Xz1 = [1 5 0]';
Yf = minreal(zpk(C*Al*Xz1));

[num_Xf,den_Xf]=tfdata(Yf,'v');
[r1,p1]=residue(num_Xf,den_Xf)
% Yt = 0

%% Punto 4
% Calcolo risposta libera (dati 2)
Xz1 = [0 0 3]';
Yf = minreal(zpk(C*Al*Xz1));

[num_Xf,den_Xf]=tfdata(Yf,'v');
[r1,p1]=residue(num_Xf,den_Xf)
% Yt = 15 * e^(-2t)

%% Punto 5
% Calcoli dei punti 2 3 4 usando ss()

% Definizione t
t = linspace(0,10,100);

% Calcolo risposta forzata con u(t) = 9 heaviside(t)
S = ss(A,B,C,0);
Yf = step(9*S,t);

% Plot risposta forzata
figure
plot(t,Yf);
title("Risposta forzata")
xlabel('t')
ylabel('Y_f')

% Calcolo risposta libera (dati 1)
Xi1 = [1 5 0]';
Y1 = initial(S, Xi1, t);

% Calcolo risposta libera (dati 2)
Xi2 = [0 0 3]';
Y2 = initial(S, Xi2, t);

% Plot due risposte
figure
tiledlayout(2,1)

nexttile
plot(t,Y1);
title("Risposta libera X(0) = [1 5 0]'")
xlabel('t')
ylabel('Y_1')

nexttile;
plot(t,Y2);
title("Risposta libera X(0) = [0 0 3]'")
xlabel('t')
ylabel('Y_2')

%% Punto bonus: somma risposta forzata con risposte libere
figure
tiledlayout(2,1)

nexttile
plot(t,Y1+Yf);
title("Risposta libera X(0) = [1 5 0]' con risposta forzata")
xlabel('t')
ylabel('Y_1f')

nexttile;
plot(t,Y2+Yf);
title("Risposta libera X(0) = [0 0 3]' con risposta forzata")
xlabel('t')
ylabel('Y_2f')