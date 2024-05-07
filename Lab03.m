%% Setup
% Inizializzazione costanti
p1 = 0.003;
p2 = 0.025;
p3 = 0.000013;
V1 = 12;
VG = 126;
n = 5/54;
Gb = 81;
Ib = 15;

% Inizializzazione solver simbolico con funzioni del sistema
syms G beta I gamma r
f = [-p1*(G-Gb)-beta*G+gamma/VG, -n*I + r/V1, -p2*beta + p3*(I-Ib)];
g = G;
x = [G, I, beta];
u = [r gamma];

%% Calcolo punto di equilibrio
% L'ingresso del sistema è dato dal testo
r_eq = 16.66667;
gamma_eq = 0;
u_eq = [r_eq gamma_eq];

%%% Soluzione a mano
% Ponendo il sistema di equazioni a zero e risolvendo a mano risulta
I_eq = r_eq/(n*V1);
beta_eq = p3/p2 .* (I_eq - Ib);
G_eq = (gamma_eq/VG + p1*Gb)/(p1+beta_eq);

x_eq = [ G_eq I_eq beta_eq ];
% Osserviamo che G_eq corrisponde al valore corretto di glicemia (81)

%%% Soluzione con solver simbolico
% Essendo il sistema non-lineare, è necessario utilizzare il motore
% simbolico per risolvere il sistema di equazioni
[G_eq, I_eq, beta_eq] = solve(subs(f,u,u_eq)==[0 0 0]);
% Il comando double calcola il valore delle frazioni
x_eq = double([ G_eq I_eq beta_eq ])

% Il risultato calcolato da MATLAB è identico a quello trovato a mano

%% Linearizzazione del sistema
% Lo jacobiano si può calcolare sia a mano che utilizzando il motore
% simbolico di MATLAB. La seguente soluzione implementa il secondo metodo

% Per ogni matrice si calcola prima lo jacobiano corrispondente e poi si
% effettua la sostituzione utilizzando il comando subs, che ci permette di
% sostituire le variabili di stato / l'input del sistema con il
% corrispondente valore nel punto di equilibrio

A = jacobian(f,x);
A = double(subs(A,x,x_eq));

B = jacobian(f,u);
B = double(subs(B,u,u_eq));

C = jacobian(g,x);
C = double(subs(C,x,x_eq));

D = jacobian(g,u);
D = double(subs(D,u,u_eq));

%% Studio stabilità
s = tf('s');
W = minreal(zpk(inv(s*eye(size(A))-A)));

%%% Stabilità interna
E = real(eig(A))
% Il sistema è stabile perché tutti gli autovalori sono negativi

%%% Stabilità BIBO
H = C*W*B;
p = pole(H)

% Il sistema + stabile BIBO in quanto tutti i poli sono negativi

%% Studio raggiungibilità e controllabilità
M_R = ctrb(A,B);
rank(M_R)
% La matrice è di rango 3, quindi il sistema è completamente raggiungibile

%% Progettazione della legge di controllo
lambda_k = [-1 -2 -3];

% Calcolo K
K = place(A,B(:,1),lambda_k);

% Calcolo alpha
alpha = inv(-(C-D(:,1)*K)*((A-B(:,1)*K)\B(:,1))+D(:,1));
