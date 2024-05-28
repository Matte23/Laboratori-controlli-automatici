%% Dati
s = tf('s');
Gp = 25/(s^3 + 3.3*s^2 + 2*s);
Gs = 1;
Ga = 0.095;
Gr = 1;
Gd = 1;

% Specifiche
Kd = 1;

%% Calcolo Kc da errore di inseguimento
% L'input di riferimento è una rampa
% Le formule sono state prese dalle tabelle sulle slide
rho_r = 0.15;
Kp = 25/2;
Kc_r = 1/(rho_r*Kp*Ga) % Kc deve essere maggiore di questo

%% Calcolo Kc da disturbo su attuatore
% L'input di riferimento è un gradino
% La formula si ottiene calcolando il limite manualmente
rho_a = 0.015;
Da0 = 0.0055;
Kc_a = Da0/(Ga*rho_a)

%% Calcolo MFLS e omega_c da disturbo dell'impianto
% Il disturbo è di tipo sinusoidale
rho_p = 0.0005;
ap = 0.02;
omega_p = 0.02;
MLFS = mag2db(rho_p/ap)
omega_l = omega_p*10.^(-MLFS/40);
omega_c_inf = 2*omega_l

%% Calcolo MHFT e omega_c da disturbo del sensore
% Il disturbo è di tipo sinusoidale
rho_s = 0.0005;
as = 0.1;
omega_s = 40;
MHFT = mag2db((rho_s*Gs)/as)
omega_h = omega_s*10.^(MHFT/40);
omega_c_inf = omega_h/2

%% Calcolo picco di risonanza e omega_c da tempo di salita e di assestamento
scap = 0.1;
tr = 3;
ts = 12;
alpha = 0.05;
xi = abs(log(scap))/sqrt(pi.^2+log(scap).^2)
Tp = 1/(2*xi*sqrt(1-xi.^2))
Sp = (2*xi*sqrt(2+4*xi^2+2*sqrt(1+8*xi^2)))/(sqrt(1+8*xi^2)+4*xi^2-1)
omega_c_rise = (((pi-acos(xi))*sqrt(sqrt(1+4*xi^4)-2*xi^2))/sqrt(1-xi^2))/tr
omega_c_settle = (-log(alpha)*sqrt(sqrt(1+4*xi^4)-2*xi^2))/(xi*ts)

% Diagramma di Nichols e Nyquist
myngridst(Tp,Sp)

%% Definizione di Gc
% Determinazione del modulo di Kc
Kc = max(Kc_a, Kc_r);

% Analisi del diagramma di Nyquist per la determinazione del segno di Kc
L = -Kc*Gp*Ga;
[num, den] = tfdata(L, "v");
figure
nyquist1(num, den)
grid on