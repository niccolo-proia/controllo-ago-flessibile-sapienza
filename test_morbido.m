% 1. Riduciamo la rigidezza del 50% (da 1000 a 500)
kt = 500 * ones(N, 1); 

% 2. Ricalcoliamo K
K = diag(6*alpha * ones(N,1) + kt);
K = K + diag(-4*alpha * ones(N-1,1), 1) + diag(-4*alpha * ones(N-1,1), -1);
K = K + diag(alpha * ones(N-2,1), 2) + diag(alpha * ones(N-2,1), -2);
K(1,1) = 7*alpha + kt(1);
K(N-1, N) = -2*alpha;
K(N-1, N-1) = 5*alpha + kt(N-1);
K(N, N-2) = 2*alpha;
K(N, N-1) = -4*alpha;
K(N, N)   = 2*alpha + kt(N);

% 3. Aggiorniamo la dinamica dell'ago (Matrice A) in memoria
A = [zeros(N, N), eye(N); 
    -inv(M)*K, -inv(M)*D];
disp('Fisica aggiornata: Tessuto MORBIDO (-50%)');