%% Phonon dispersion of a 1D diatomic chain with alternating spring 
% constants C, 10C

%% Physical parameters
k1 = 1; % spring 1
k2 = 2; % spring 2
m1 = 5; % mass 1
m2 = 1; % mass 2
a  = 1; % lattice constant 

% k-space 
N = 100;
k = linspace(-pi/a, pi/a, N);

% Initialize to store angular frequencies
w_a = zeros(1, N);
w_o  = zeros(1, N);

for n = 1:N
    ki = k(n);

    %% Matrix Components of A (System is Au = w^2Mu)
    A1 = k1 + k2;
    A2 = -(k1 + k2*exp(-1i*ki*a));
    A3 = conj(A2);
    A4 = k1 + k2;

    % Matrix A
    A = [A1,A2;
         A3,A4];

    % Mass matrix
    M = [m1,0;
         0,m2];

    % Solving for w^2 (Du = ω^2Mu)
    w2 = eig(M \ A);  
    w2 = sort(real(w2));

    % Angular frequencies (higher freq is optical band, lower freq is acou)
    w_a(n) = sqrt(w2(1));
    w_o(n)  = sqrt(w2(2));
end

%% Plotting dispersion
clf; 
plot(k, w_a, 'k.'); hold on;
plot(k, w_o, 'r.');
xlabel('Wavevector k (rad/m)');
ylabel('Angular Frequency \omega (rad/s)');
title('Phonon Dispersion of Diatomic Basis');
legend('Acoustic branch', 'Optical branch');

