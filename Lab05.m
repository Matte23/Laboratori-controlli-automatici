%% Setup variabili
A = [0 1; -3 -4];
B = [0 ; 1];
C = [2 1];
D = 0;

%% Punto 1-2
lambda_k = [-0.7 -0.8];
lambda_o = [-10 -11];

% Controllo raggiungibilità
M_R = ctrb(A,B);
rank(M_R)
% Rango 2: completamente raggiungibile

% Controllo osservabilità
M_O = obsv(A,C);
rank(M_O)
% Rango 2: completamente osservabile

% Calcolo K
K = place(A,B,lambda_k);

% Calcolo L
L = place(A',C',lambda_o)';

% Calcolo alpha
alpha = inv(-C*((A-B*K)\B));

%% Punto 3
% (Su simulink) Si osserva che l'uscita del sistema tende ad 1

%% Punto 4
lambda_k_4 = [-10 -12];
K_4 = place(A,B,lambda_k_4);
alpha_4 = inv(-C*((A-B*K_4)\B));
% Si osserva che con autovalori troppo elevati il sistema converge lo stesso, ma
% all'inizio si osserva un transitorio con valori troppo grandi (più sono grandi gli autovalori, più
% l'errore cresce)

% Con questa tecnica non riesco a controllare il sistema velocemente per
% via dello zero nella H(s)
minreal(zpk(ss(A,B,C,D)))

%% Punto 5
% Nota: ripristinare gli autovalori del punto 1.
epsilon = 10.^(-2/20);
B_reale = B.*epsilon;
% Osserviamo che il sistema "reale" (con la perturbazione) è instabile.

%% Punto 6
lambda_o_6 = [-1000 -1100];
L_6 = place(A',C',lambda_o_6)';

%% Punto 7
lambda_o_7 = [-1000 -2];
L_7 = place(A',C',lambda_o_7)';