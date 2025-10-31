%% Periodic Crystal Lattice in 1D
% Periodic boundary conditions, u_n+N = u_n, chain of atoms of length Na,
% approximating the interactions as harmonic, springs connecting atoms

%% Parameters (for a simple cubic crystal of Polonium)
N = 100; % Number of atoms
a = 3e-10; % Lattice constant in meters (3 Angstroms here)
m = (209)*1.66e-27; % Atomic mass (kg), (where 1.66e-27 is kg per amu)
k_spring = 100; % Spring constant between atoms (N/m)

% Constants
k_B = 1.38e-23; % Stefan-Boltzmann constant (J/K)
hbar = 1.05e-34;  % Planck constant over 2pi, h/2pi (J*s)

%% Wavevectors (k's) and Corresponding Frequencies (w's)
% (Only those in the 1st Brullion Zone (B.Z.) are needed 
% because these are the only independent wavevectors, the reciprocal 
% lattice is periodic in reciprocal lattice vectors, g = 2pi/a)

n = -N/2:(N/2); % We only need wavevectors in the 1st B.Z. (from -pi/a to 
% pi/a), these values of n give the correct values of k below
k = 2*pi*n/(N*a); % Periodic B.C., so u_n+N = u_n, look for solutions of 
% form u_n = Aexp(i(kx_n-wt)), valid for small oscillations, where the 
% position of the nth atom is x_n = na. This gives the condition on k
% where exp(ikNa) = 1 so k = 2pi*n/Na. 

w = sqrt(4*k_spring/m)*abs(sin(k*a/2)); % Phonon dispersion equation, 
% obtained by solving for the equation of motion for a 1D lattice of atoms
% connected by springs, with a spacing between them of a (the lattice
% constant).

%% Plotting Dispersion (w vs k)
figure(1); clf;
plot(k, w, '.k');
xlabel('k (m^{-1})');
ylabel('\omega (radians/s)');
title('Phonon Dispersion in 1D');


%% Heat Capacity
T = linspace(1, 1000, 100);  % Temperature in Kelvin
Cv = zeros(size(T)); % Preallocate
for i = 1:length(T)
    x = hbar*w./(k_B*T(i)); % Factor in exponent of specific heat
    x(w==0) = 0.001; % To not allow division by 0 
    Cv(i) = k_B * sum( (x.^2 .* exp(x)) ./ (exp(x) - 1).^2 ) / length(w);
    % Heat capacity formula using the energy of quantum harmonic
    % oscillators E_n = hbar*w(n + 1/2) for the phonons and thermodynamic
    % partition function Z = sum_n(exp(-betaE_n)) where beta is 1/(k_B*T) 
    % which can be used to find the average energy <E> = -dln(Z)/d(beta),
    % which can be used to find specific heat at constant volume with 
    % C_v = d<E>/dT = -(1/(k_B*T^2))(d<E>/d(beta)). Also this sum is over
    % all vibrational modes, all w's, so the number of modes (w's) is
    % divided out to normalize.
end
Cv_norm = Cv/k_B; % Normalize to High Temperature (Dulong-Petit) limit 
% (1 factor of k_B in 1D), in 3D there would be 3 factors (in classical
% calculation, equipartition theorem gives (1/2)k_B*T of energy for each
% atom and 3 degrees of freedom in 3D so internal energy is (3/2)*N*k_B*T,
% for N atoms, so the heat capacity is (3/2)*N*k_B)

% Plotting Phonon contribution to heat capacity in 1D
figure(2); clf;
plot(T, Cv_norm, '.r');
xlabel('Temperature (K)');
ylabel('Normalized Heat Capacity C_v/k_B');
title('Phonon Heat Capacity in 1D');

