%% Setup
% Inizializzazione costanti
m = 0.02;
G = 9.81;
Kt = 708.27;
Km = 1.52*10^-4;

% Inizializzazione solver simbolico con funzioni del sistema
syms x1 x2 Im
f = [x2, G-(Km.*Im.^2)/(m.*x1.^2)];
x = [x1, x2];
u = Im;

%% Calcolo punto di equilibrio
% L'ingresso del sistema è dato dal testo
u_eq = 0.8;

% La soluzione si può sia calcolare a mano che usando il solver simbolico
% di matlab
[x1_eq, x2_eq] = solve(subs(f,u,u_eq)==[0 0]);
% Il comando double calcola il valore delle frazioni
x_eq = double([ x1_eq, x2_eq ]);

% Si osserva che sono presenti due punti di equilibrio. Prendiamo in
% considerazione solo quello con variabili positive
x_eq = x_eq(2,:)

%% Linearizzazione del sistema
% Lo jacobiano si può calcolare sia a mano che utilizzando il motore
% simbolico di MATLAB. La seguente soluzione implementa il secondo metodo

% Definisco la funzione uscita attorno al punto di equilibrio
g = Kt.*(x1-x_eq(1));

A = jacobian(f,x);
A = double(subs(A,[x u],[x_eq u_eq]));

B = jacobian(f,u);
B = double(subs(B,[x u],[x_eq u_eq]));

C = jacobian(g,x);
C = double(subs(C,x,x_eq));

D = jacobian(g,u);
D = double(subs(D,u,u_eq));

%% Studio stabilità
s = tf('s');
W = minreal(zpk(inv(s*eye(size(A))-A)));

%%% Stabilità interna
E = real(eig(A))
% Il sistema è instabile perché è presente un autovalore positivo

%%% Stabilità BIBO
H = C*W*B;
p = pole(H)

% Il sistema è instabile BIBO è presente un polo positivo

%% Studio raggiungibilità e osservabilità
% Raggiungibilità
M_R = ctrb(A,B);
rank(M_R)
% Il rango è due, quindi il sistema è completamente raggiungibile

% Osservabilità
M_O = obsv(A,C);
rank(M_O)
% Il rango è due, quindi il sistema è completamente osservabile

%% Progettazione sistema di controllo
% Scelta autovalori
lambda_k = [-0.7 -0.8];
lambda_o = [-10 -11];

% Calcolo K
K = place(A,B,lambda_k);

% Calcolo L
L = place(A',C',lambda_o)';

% Calcolo alpha
alpha = inv(-C*((A-B*K)\B));
