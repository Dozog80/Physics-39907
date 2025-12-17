%% Periodic Potential Band Structure
% Physical constants
hbar = 1;
m = 1;
a = 1;
G = 2*pi/a;

% Fourier coefficients of the periodic potential
% V(x) = sum_G (V_G*exp(iGx))
% including ±G, ±2G, ±3G
Vcoeffs = struct();
Vcoeffs(1).G = 1*G;   Vcoeffs(1).V = 7;
% Vcoeffs(2).G = 2*G;   Vcoeffs(2).V = 1;
% Vcoeffs(3).G = 3*G;   Vcoeffs(3).V = 1/2;

% Number of Fourier components (plane waves)
M = max([Vcoeffs.G] / G); % max harmonic index
Npw = 2*M + 1; % -M to M

% k-space inside 1st BZ
Nk = 100;
k = linspace(-G/2, G/2, Nk);

% Free-electron energy
E0 = @(q) (hbar^2 * q.^2)/(2*m);
Ebands = zeros(Npw, Nk); % Initialize

%% Build Hamiltonian for each k
for i = 1:Nk
    ki = k(i);

    % reciprocal lattice shifts, -MG to MG
    m_list = (-M:M);
    kstates = ki + m_list*G;
    % Initialize Hamiltonian
    H = zeros(Npw, Npw);

    % Diagonal terms (free electron energies)
    for p = 1:Npw
        H(p,p) = E0(kstates(p));
    end

    % Off-diagonal terms from periodic potential
    for c = 1:length(Vcoeffs)
        Gc = Vcoeffs(c).G;
        Vc = Vcoeffs(c).V;

        for p = 1:Npw
            for q = 1:Npw
                if abs(kstates(p) - kstates(q) - Gc) < 1e-12
                    H(p,q) = H(p,q) + Vc;
                end
                if abs(kstates(p) - kstates(q) + Gc) < 1e-12
                    H(p,q) = H(p,q) + conj(Vc);
                end
            end
        end
    end

    % Diagonalize Hamiltonian
    Evals = sort(eig(H));
    Ebands(:,i) = Evals;
end

% Fermi energy, 1D system
N_e = 50; % number of electrons
Lbox = 100; % length
n = N_e / Lbox;
kF = pi*n;
EF = E0(kF);

%% Plot bands (Reduced Zone Scheme)
figure(1); clf; hold on; box on;
for band = 1:Npw
    plot(k, Ebands(band,:), 'k.');
end

% Fermi energy
yline(EF, 'k--');
text(0, EF + 3, sprintf('E_F = %.3f', EF));
xlabel('k (rad/m)'); ylabel('E (Energy eV)');
title(sprintf('Electron in Periodic Potential Band Structure \n (Reduced Zone Scheme)'));
