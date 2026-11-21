% 1. Pesi unitari per filtrare il rumore ad alta frequenza
Q_aug = eye(41);

% Ripristiniamo l'integratore per dare al sistema la forza di raggiungere il bersaglio
Q_aug(41, 41) = 10; 

% Il compromesso d'oro: abbastanza alto da smorzare l'overshoot, 
% abbastanza basso da permettere il movimento entro i 5 secondi.
R_aug = 50; 

% 2. Sintesi LQR CONTINUA
K_aug = lqr(A_aug, B_aug, Q_aug, R_aug);

% 3. Separazione per Simulink
K_x = K_aug(1:40);
K_i = K_aug(41);