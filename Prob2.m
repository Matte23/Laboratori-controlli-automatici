%% Dati
s = tf('s');
Gp = 40/(s^2 + 3*s + 4.5);
Gs = 1;
Ga = -0.09;
Gr = 1;
Gd = 1;

% Specifiche
Kd = 1;
p = 0; % L'impianto non ha poli nell'origine

%% Calcolo Kc da errore di inseguimento
% L'input di riferimento è una rampa
% Le formule sono state prese dalle tabelle sulle slide

h_r = 1; % input è una rampa
nu_r = h_r - p % limite inferiore per nu_r

rho_r = 0.35;
Kp = evalfr(Gp, 0);
Kc_r = -1/(rho_r*Kp*Ga) % Kc deve essere maggiore di questo

%% Calcolo Kc da disturbo su attuatore
% Il disturbo è un gradino
h_a = 0;
nu_a = h_a - p
% da questo calcolo risulta nu>=0, ma considerando che il punto precedente 
% impone nu >= 1, allora risulta calcolando il limite che il contributo del
% disturbo tende a zero, pertanto non possiamo ricavare nessuna condizione
% su Kc da questo requisito

%% Calcolo Kc da disturbo dell'impianto
% Il disturbo è una rampa
rho_p = 0.001;
Dp0 = 0.003;
Kc_p = -Dp0/(rho_p*Kp*Ga)

%% Calcolo MHFT e omega_c da disturbo del sensore
% Il disturbo è di tipo sinusoidale
rho_s = 0.0002;
as = 0.01;
omega_s = 50;
MHFT = mag2db((rho_s*Gs)/as)
omega_h = omega_s*10.^(MHFT/40);
omega_c_inf = omega_h/2

%% Calcolo picco di risonanza e omega_c da tempo di salita e di assestamento
scap = 0.08;
tr = 2.5;
ts = 10;
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
Kc = max(Kc_a, Kc_p);

% Analisi del diagramma di Nyquist per la determinazione del segno di Kc
L = -Kc*Gp*Ga;
[num, den] = tfdata(L, "v");
figure
nyquist1(num, den)
grid on