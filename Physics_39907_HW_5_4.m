N = 20; % Number of springs, 
c = 1 + (1:N)/2; % Spring constants c_m = 1 + m/2
C = diag(c);  % Diagonal spring constants matrix
M = N - (1:N-1); % Masses
f = M'; % Forces

Delta = toeplitz([1, -1, zeros(1, N-2)], [1, zeros(1, N-2)]); % The First
% Differences operator

K = Delta' * C * Delta; % The definition of K
u = K \ f % Solving the system, shows displacements for each mass
residual = K*u - f % To see floating point errors

% Plotting Displacements versus n
clf;
n = 1:N-1;
plot(n, u, 'o-b'); % Plot
xlabel('n'); % Index n
ylabel('Displacements u_n'); % Displacements
title('Displacements u_n for N = 20');

