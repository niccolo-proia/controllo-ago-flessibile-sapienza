% =========================================================================
% Tesi: Modellazione di un ago flessibile (Trave di Eulero-Bernoulli)
% Metodo: Differenze Finite (FDM)
% =========================================================================

clear all; close all; clc;

%% 1. Definizione dei Parametri Fisici (Ecco qui i valori reali)
L = 0.15;           % Lunghezza dell'ago [m]
N = 20;             % Numero di nodi per la discretizzazione spaziale
h = L / N;          % Passo di discretizzazione spaziale [m]

E_young = 200e9;    % Modulo di Young dell'ago (es. acciaio) [Pa]
I_mom = 1e-12;      % Momento d'inerzia della sezione trasversale [m^4]
rho = 8000;         % Densità del materiale dell'ago [kg/m^3]
Area = 1e-6;        % Area della sezione trasversale [m^2]
c = 0.5;            % Coefficiente di smorzamento viscoso [N*s/m]

% Parametri derivati
alpha = (E_young * I_mom) / (h^4);
m = rho * Area;     % Massa per unità di lunghezza

% Profilo della fondazione elastica alla Winkler (Rigidezza del tessuto)
% Per ora lo considero costante, ma può essere una funzione di x
kt = 1000 * ones(N, 1); 

%% 2. Costruzione delle Matrici Massa, Smorzamento, Rigidezza (M, D, K)
M = m * eye(N);     % Matrice di Massa (diagonale)
D = c * eye(N);     % Matrice di Smorzamento (diagonale)

% Costruzione della Matrice di Rigidezza K (matrice bandata)
% Inizializzo la diagonale principale (6*alpha + kt)
K = diag(6*alpha * ones(N,1) + kt);

% Aggiungo le prime codiagonali (-4*alpha)
K = K + diag(-4*alpha * ones(N-1,1), 1) + diag(-4*alpha * ones(N-1,1), -1);

% Aggiungo le seconde codiagonali (alpha)
K = K + diag(alpha * ones(N-2,1), 2) + diag(alpha * ones(N-2,1), -2);

% APPLICAZIONE DELLE CONDIZIONI AL CONTORNO
% Nodo 1 (Incastro mobile, u(t))
K(1,1) = 7*alpha + kt(1);

% Nodo N (Estremo libero, forza di bisello Fp)
K(N-1, N) = -2*alpha;
K(N-1, N-1) = 5*alpha + kt(N-1);
K(N, N-2) = 2*alpha;
K(N, N-1) = -4*alpha;
K(N, N)   = 2*alpha + kt(N);

%% 3. Costruzione dei Vettori degli Ingressi (Bu, Ed)
Bu = zeros(N, 1);
Bu(1) = 4*alpha;
Bu(2) = -alpha;

Ed = zeros(N, 1);
Ed(N) = 2/h;

%% 4. Costruzione del Modello in Spazio di Stato (A, B, C, E)
% Matrici di zeri e identità per facilitare la costruzione a blocchi
Z = zeros(N, N);
I_mat = eye(N);

% Matrice A (Dinamica del sistema)
A = [Z, I_mat;
    -inv(M)*K, -inv(M)*D];

% Matrici B ed E (Ingressi sulle accelerazioni)
B = [zeros(N, 1);
    inv(M)*Bu];

E_dist = [zeros(N, 1);
    inv(M)*Ed];

% Matrice C (Misura della deflessione di punta)
C = zeros(1, 2*N);
C(N) = 1;           % L'uscita è il nodo N (la punta)

%% 5. Creazione del Sistema LTI in MATLAB
% Il sistema ha 2 ingressi: u(t) alla base e il disturbo Fp(t) in punta.
% Affianco B ed E per creare la matrice d'ingresso completa [B E_dist]
% La matrice D (feedthrough) è nulla, quindi mettiamo [0 0]
sys = ss(A, [B, E_dist], C, [0, 0]);

% Rinomino gli ingressi e le uscite per chiarezza nei grafici
sys.InputName = {'Spostamento base u(t)', 'Forza bisello F_p'};
sys.OutputName = {'Deflessione di punta w(L,t)'};

disp('Modello in spazio di stato creato con successo!');

% Tempo di campionamento dell'ecografo (es. 30 Hz -> 0.033 secondi)
Ts = 0.033; 

% Converte il sistema continuo in discreto mantenendo la precisione
sys_d = c2d(sys, Ts, 'zoh');

% Mappa Poli-Zeri in anello aperto (L = 0.15 m)
figure;
pzmap(sys);
grid on;
title('Mappa Poli-Zeri del Sistema Ago-Tessuto (L = 0.15 m)');

%% 6. Analisi delle Forme Modali
% Calcolo degli autovalori (frequenze quadrate) e autovettori (forme modali)
[V, Omega2] = eig(K, M);

% Estrazione delle frequenze (ordinandole dalla più bassa alla più alta)
[omega2_ordinate, indici] = sort(diag(Omega2));
V_ordinati = V(:, indici);

% Plot delle prime 3 forme modali
figure;
hold on;
asse_x = linspace(0, L, N);
plot(asse_x, V_ordinati(:,1) / max(abs(V_ordinati(:,1))), 'r', 'LineWidth', 2);
plot(asse_x, V_ordinati(:,2) / max(abs(V_ordinati(:,2))), 'g', 'LineWidth', 2);
plot(asse_x, V_ordinati(:,3) / max(abs(V_ordinati(:,3))), 'b', 'LineWidth', 2);
title('Prime 3 Forme Modali dell''Ago Flessibile');
xlabel('Posizione lungo l''ago (m)');
ylabel('Deformazione Normalizzata');
legend('Modo 1', 'Modo 2', 'Modo 3');
grid on;

%% 6b. Effetto della Profondità di Inserzione (L = 0.075 m) - Variabili isolate
L_ridotto = 0.075;
h_ridotto = L_ridotto / N;
alpha_ridotto = (E_young * I_mom) / (h_ridotto^4);

K_ridotto = diag(6*alpha_ridotto * ones(N,1) + kt) ...
    + diag(-4*alpha_ridotto * ones(N-1,1), 1) + diag(-4*alpha_ridotto * ones(N-1,1), -1) ...
    + diag(alpha_ridotto * ones(N-2,1), 2) + diag(alpha_ridotto * ones(N-2,1), -2);
K_ridotto(1,1) = 7*alpha_ridotto + kt(1);
K_ridotto(N-1, N) = -2*alpha_ridotto;
K_ridotto(N-1, N-1) = 5*alpha_ridotto + kt(N-1);
K_ridotto(N, N-2) = 2*alpha_ridotto;
K_ridotto(N, N-1) = -4*alpha_ridotto;
K_ridotto(N, N)   = 2*alpha_ridotto + kt(N);

A_ridotto = [Z, I_mat; -inv(M)*K_ridotto, -inv(M)*D];
sys_ridotto = ss(A_ridotto, [B, E_dist], C, [0, 0]);

figure;
pzmap(sys_ridotto);
grid on;
title('Mappa Poli-Zeri a Profondità Dimezzata (L = 0.075 m)');
%% 7. Sintesi del Controllore PID
% Isolo la dinamica ingresso-punta (dal primo ingresso u(t) all'unica uscita y)
sys_u = sys(1, 1); 

% Progetto il controllore PID in automatico
% pidtune cerca il miglior compromesso tra tempo di assestamento e robustezza
[C_pid, info] = pidtune(sys_u, 'PID');

% Mostro a video i parametri calcolati e le prestazioni teoriche
disp('Parametri e prestazioni del controllore PID:');
C_pid
info

% Costruisco il sistema in anello chiuso per l'inseguimento del riferimento
% feedback(Azione_Diretta, Azione_Retroazione)
sys_cl_pid = feedback(sys_u * C_pid, 1);

% Test: Risposta al gradino del sistema ad anello chiuso (riferimento = 1 mm)
figure;
step(sys_cl_pid);
title('Inseguimento del bersaglio con Controllore PID');
ylabel('Deflessione di punta (mm)');
grid on;

%% 7b. Reiezione del Disturbo e Sforzo di Controllo (PID)

% 1. REIEZIONE DEL DISTURBO (Forza di bisello F_p)
% Estraggo la dinamica dal disturbo (ingresso 2) all'uscita (punta)
sys_dist = sys(1, 2); 

% Calcolo la funzione di trasferimento ad anello chiuso per il disturbo
% Formula teorica: Y(s)/D(s) = G_d(s) / (1 + C(s)*G_u(s))
sys_cl_dist = sys_dist / (1 + C_pid * sys_u);

% Test: Risposta al gradino del disturbo (applicazione di 1 N in punta)
figure;
step(sys_cl_dist);
title('Reiezione del disturbo di bisello (F_p) con Controllore PID');
ylabel('Deflessione di punta (m)');
grid on;

% 2. SFORZO DI CONTROLLO (Azione sulla base u(t))
% Funzione di trasferimento dal riferimento all'ingresso di controllo u
% Formula teorica: U(s)/R(s) = C(s) / (1 + C(s)*G_u(s))
sys_cl_u = feedback(C_pid, sys_u);

figure;
step(sys_cl_u);
title('Sforzo di controllo del PID (Spostamento alla base)');
ylabel('Spostamento u(t) (mm)');
grid on;


%% 8. Sintesi del Controllore Ottimo LQR con Integratore

% 1. CREAZIONE DEL MODELLO AUMENTATO
% Isolo la matrice B relativa al solo attuatore alla base (prima colonna)
Bu_base = sys.B(:, 1);

% Dimensioni correnti dello stato (2*N = 40)
n_stati = size(sys.A, 1);

% Costruisco le matrici aumentate (41x41 e 41x1)
A_aug = [sys.A, zeros(n_stati, 1); 
    sys.C, 0];

B_aug = [Bu_base; 
    0];

% 2. DEFINIZIONE DELLE MATRICI DI PESO (Q e R)
% LQR minimizza il costo J = Somma(x'*Q*x + u'*R*u)
% Q pesa gli stati (più è alto, più il sistema cerca di annullare l'errore rapidamente)
% R pesa lo sforzo di controllo (più è alto, più il sistema risparmia energia)

Q = eye(n_stati + 1);      % Inizializzo Q come matrice identità 41x41
Q(end, end) = 1e6;         % Penalizzo enormemente l'errore integrale (stato 41) per forzare l'inseguimento esatto

R = 1;                     % Costo unitario per lo sforzo di controllo u(t)

% 3. CALCOLO DEI GUADAGNI OTTIMI
% La funzione lqr calcola il vettore dei guadagni K_aug (1x41)
K_aug = lqr(A_aug, B_aug, Q, R);

% Separo i guadagni: 
% I primi 40 moltiplicheranno lo stato stimato dall'osservatore (K_x)
% L'ultimo moltiplicherà l'errore integrale (K_i)
K_x = K_aug(1:n_stati);
K_i = K_aug(end);

disp('Guadagni LQR calcolati con successo (K_x e K_i).');


%% 8b. Simulazione dell'Anello Chiuso LQR e Confronto

% 1. COSTRUZIONE DEL SISTEMA AD ANELLO CHIUSO LQR
% Le nuove equazioni della dinamica chiusa con integratore (u = -K_x*x - K_i*x_i)
A_cl_lqr = [sys.A - Bu_base*K_x, -Bu_base*K_i; 
    sys.C,               0];

% L'ingresso del sistema chiuso è il riferimento r(t)
B_cl_lqr = [zeros(n_stati, 1); 
    -1];

% L'uscita rimane la deflessione in punta (primi 40 stati)
C_cl_lqr = [sys.C, 0];

D_cl_lqr = 0;

% Creo l'oggetto di stato per l'LQR chiuso
sys_cl_lqr = ss(A_cl_lqr, B_cl_lqr, C_cl_lqr, D_cl_lqr);

% 2. TEST E GRAFICI: INSEGUIMENTO DEL RIFERIMENTO
figure;
step(sys_cl_lqr);
title('Inseguimento del bersaglio con Controllore LQR (Integratore)');
ylabel('Deflessione di punta (mm)');
grid on;

% 3. SFORZO DI CONTROLLO LQR
% Calcolo l'azione u(t) generata dall'LQR per raggiungere il bersaglio
C_u_lqr = [-K_x, -K_i];
sys_u_lqr = ss(A_cl_lqr, B_cl_lqr, C_u_lqr, 0);

figure;
step(sys_u_lqr);
title('Sforzo di controllo LQR (Spostamento alla base)');
ylabel('Spostamento u(t) (mm)');
grid on;


%% 8c. Analisi di Sensibilità dei Pesi LQR (Q e R)

% Preparo le tre configurazioni di pesi (Q_int e R)
% 1. Standard (Quello che hai appena calcolato)
Q_int_1 = 1e6;   R_1 = 1;
% 2. Aggressivo (Penalizzo moltissimo l'errore, non mi importa dell'energia)
Q_int_2 = 1e8;   R_2 = 1;
% 3. Lento/Risparmioso (Penalizzo molto l'uso del motore, accetto un ritardo)
Q_int_3 = 1e6;   R_3 = 100;

% Vettori per salvare i sistemi chiusi e tracciarli insieme
sys_cl_list = cell(1,3);
sys_u_list = cell(1,3);
config_names = {'LQR Standard', 'LQR Aggressivo', 'LQR Risparmio Energia'};

pesi_Q = [Q_int_1, Q_int_2, Q_int_3];
pesi_R = [R_1, R_2, R_3];

for i = 1:3
    % Aggiorno i pesi
    Q_temp = eye(n_stati + 1);
    Q_temp(end, end) = pesi_Q(i);
    R_temp = pesi_R(i);

    % Calcolo i nuovi guadagni
    K_aug_temp = lqr(A_aug, B_aug, Q_temp, R_temp);
    K_x_temp = K_aug_temp(1:n_stati);
    K_i_temp = K_aug_temp(end);

    % Costruisco l'anello chiuso
    A_cl_temp = [sys.A - Bu_base*K_x_temp, -Bu_base*K_i_temp; 
        sys.C,                    0];

    % Salvo la funzione di trasferimento per l'inseguimento (punta)
    sys_cl_list{i} = ss(A_cl_temp, B_cl_lqr, C_cl_lqr, 0);

    % Salvo la funzione di trasferimento per lo sforzo di controllo (base)
    sys_u_list{i} = ss(A_cl_temp, B_cl_lqr, [-K_x_temp, -K_i_temp], 0);
end

% Genero i grafici a confronto
figure;
step(sys_cl_list{1}, sys_cl_list{2}, sys_cl_list{3});
title('Confronto Inseguimento Bersaglio al variare di Q e R');
ylabel('Deflessione di punta (mm)');
legend(config_names);
grid on;

figure;
step(sys_u_list{1}, sys_u_list{2}, sys_u_list{3});
title('Confronto Sforzo di Controllo al variare di Q e R');
ylabel('Spostamento u(t) (mm)');
legend(config_names);
grid on;


%% 9. Sintesi dell'Osservatore di Luenberger (Obiettivo O4)

% Il sistema ha 40 stati (posizioni e velocità dei 20 nodi)
% Per evitare errori numerici di malcondizionamento dovuti alla stiffness,
% usiamo la dualità LQR per calcolare la matrice dei guadagni dell'osservatore (L_obs).

% Q_obs pesa la rapidità con cui vogliamo che l'osservatore corregga l'errore di stima
Q_obs = eye(n_stati) * 1e12; 

% R_obs pesa l'affidabilità del sensore (l'ecografo)
R_obs = 1;

% Calcolo dei guadagni dell'osservatore (L_obs ha dimensione 40x1)
L_obs = lqr(sys.A', sys.C', Q_obs, R_obs)';

% Costruzione della dinamica dell'osservatore
% Equazione: dx_est = A*x_est + B*u + L_obs*(y - C*x_est)
% Raggruppando per x_est: dx_est = (A - L_obs*C)*x_est + [B_u, L_obs] * [u; y]

A_obs = sys.A - L_obs * sys.C;

% L'osservatore prende in ingresso il comando del motore (u) e la misura letta (y)
B_obs = [Bu_base, L_obs];

% L'uscita dell'osservatore è l'intero vettore di stato stimato (tutti i 40 stati)
C_obs = eye(n_stati); 
D_obs = zeros(n_stati, 2);

% Genero l'oggetto spazio di stato per il solo Osservatore
sys_luenberger = ss(A_obs, B_obs, C_obs, D_obs);

disp('Guadagni dell''Osservatore di Luenberger calcolati con successo.');


%% 9b. Simulazione e Convergenza dell'Osservatore di Luenberger

% 1. CREAZIONE DELLE CONDIZIONI INIZIALI
% x0_reale: L'ago fisico è già flesso (es. deflessione lineare fino a 5 cm in punta)
x0_reale = zeros(n_stati, 1);
x0_reale(1:N) = linspace(0, 0.05, N)'; 

% x0_stimato: L'osservatore parte "cieco", credendo che l'ago sia fermo e dritto a zero
x0_stimato = zeros(n_stati, 1);       

% 2. COSTRUZIONE DEL SISTEMA COMBINATO (Realtà + Stima)
% Unisco la dinamica dell'ago vero (sys.A) e quella dell'osservatore (A_obs)
% Il termine L_obs*sys.C rappresenta l'iniezione dell'errore di misura dalla punta reale all'osservatore
A_tot = [sys.A,        zeros(n_stati);
    L_obs*sys.C,  A_obs];

B_tot = zeros(2*n_stati, 1); % Simuliamo un'evoluzione libera (nessun motore u(t) attivo)

% 3. ESTRAZIONE DI UN NODO NASCOSTO (es. il nodo centrale dell'ago)
nodo_centrale = round(N/2); % Per N=20, questo è il nodo 10

% Voglio vedere sul grafico solo 2 curve: 
% Riga 1: la posizione REALE del nodo centrale
% Riga 2: la posizione STIMATA del nodo centrale
C_tot = zeros(2, 2*n_stati);
C_tot(1, nodo_centrale) = 1;              % Estrae w_10 reale (tra i primi 40 stati)
C_tot(2, n_stati + nodo_centrale) = 1;    % Estrae w_10 stimato (tra i secondi 40 stati)

D_tot = zeros(2, 1);

% Oggetto di stato per la simulazione
sys_tot = ss(A_tot, B_tot, C_tot, D_tot);

% 4. TEST DI CONVERGENZA (Grafico Sovrapposto corretto)
X0_tot = [x0_reale; x0_stimato]; % Vettore iniziale da 80 elementi

% Definisco un vettore temporale per la simulazione (da 0 a 0.3 secondi)
t_sim = 0:0.001:0.3;

% Simulo ed estraggo i dati (Y conterrà due colonne: Out1 e Out2)
[Y_sim, t_sim, X_sim] = initial(sys_tot, X0_tot, t_sim);

figure;
plot(t_sim, Y_sim(:, 1), 'b', 'LineWidth', 1.5); % Posizione Reale (Blu)
hold on;
plot(t_sim, Y_sim(:, 2), 'r--', 'LineWidth', 1.5); % Posizione Stimata (Rossa tratteggiata)
title(sprintf('Convergenza Luenberger: Nodo nascosto %d', nodo_centrale));
xlabel('Tempo (secondi)');
ylabel('Deflessione (m)');
legend('Posizione Reale (Invisibile all''ecografo)', 'Posizione Stimata (Osservatore)');
grid on;


%% 10. Sintesi del Filtro di Kalman

% Matrice G: indica come il rumore di processo entra negli stati
G_k = eye(n_stati);

% Q_k: Covarianza del rumore di processo (incertezza del modello/tessuto)
Q_k = eye(n_stati) * 1e-1; 

% R_k: Covarianza del rumore dell'ecografo 
% (es. deviazione standard di 1 mm -> varianza di 1e-6)
R_k = 1e-6; 

% Sintesi del guadagno ottimo di Kalman (L_k)
L_k = lqe(sys.A, G_k, sys.C, Q_k, R_k);

% Costruzione della dinamica del Filtro di Kalman
A_kalman = sys.A - L_k * sys.C;
B_kalman = [Bu_base, L_k];
C_kalman = eye(n_stati);
D_kalman = zeros(n_stati, 2);
sys_kalman = ss(A_kalman, B_kalman, C_kalman, D_kalman);

disp('Filtro di Kalman calcolato con successo.');

%% 10b. Sfida Finale: Luenberger vs Kalman con Rumore Ecografico

% Vettore temporale
t_kal = 0:0.001:0.3;

% 1. LA REALTA' FISICA (L'ago vibra partendo da 5 cm)
[Y_vera_punta, ~, X_veri] = initial(sys, x0_reale, t_kal); 
Pos_Reale_Nodo10 = X_veri(:, nodo_centrale);

% 2. IL SENSORE (Aggiungo rumore bianco Gaussiano alla punta)
rumore_sensore = sqrt(R_k) * randn(length(t_kal), 1);
Y_misurata_punta = Y_vera_punta + rumore_sensore;

% 3. LA GARA DEGLI ESTIMATORI (Entrambi partono ciechi da zero)
% Ingressi: Motore fermo (0) e misura ecografica rumorosa
ingressi_stim = [zeros(length(t_kal), 1), Y_misurata_punta];

% Luenberger subisce il segnale rumoroso
[X_luenberger, ~] = lsim(sys_luenberger, ingressi_stim, t_kal, x0_stimato);
Pos_Luenberger_Nodo10 = X_luenberger(:, nodo_centrale);

% Kalman filtra il segnale rumoroso
[X_kalman, ~] = lsim(sys_kalman, ingressi_stim, t_kal, x0_stimato);
Pos_Kalman_Nodo10 = X_kalman(:, nodo_centrale);

% 4. GRAFICO DI CONFRONTO
figure;
plot(t_kal, Pos_Reale_Nodo10, 'b', 'LineWidth', 1.5);
hold on;
plot(t_kal, Pos_Luenberger_Nodo10, 'r', 'LineWidth', 0.5);
plot(t_kal, Pos_Kalman_Nodo10, 'g', 'LineWidth', 2);
title(sprintf('Ricostruzione del Nodo %d con Rumore Ecografico', nodo_centrale));
xlabel('Tempo (secondi)');
ylabel('Deflessione (m)');
legend('Realtà (Invisibile)', 'Luenberger (Vittima del rumore)', 'Kalman (Filtrato)');
grid on;