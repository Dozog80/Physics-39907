%% 3D Simple Cubic Lattice Fermi Surface
a = 1; % lattice constant
G0 = 2*pi/a;% Magnitude of Reciprocal Lattice Vector

% Using 6 reciprocal lattice vectors ±Gx, ±Gy, ±Gz
Gvecs = [
    1  0  0;
    -1  0  0;
    0  1  0;
    0 -1  0;
    0  0  1;
    0  0 -1
    ] * G0;

% Fourier coefficients of potential
Vmag = 4.5; % potential amplitude
V_G  = Vmag * ones(size(Gvecs,1),1);
NPW = size(Gvecs,1) + 1; % Plane waves, including G=0


%% 3D k-space (1st BZ)
Nk = 40; % Number of points on grid
k1d = linspace(-pi/a, pi/a, Nk);
[kx, ky, kz] = ndgrid(k1d, k1d, k1d);


%% Bands
hbar = 1;  m = 1;
Ebands = zeros(Nk, Nk, Nk); % Initialize

for ix = 1:Nk
    for iy = 1:Nk
        for iz = 1:Nk
            k = [kx(ix,iy,iz), ky(ix,iy,iz), kz(ix,iy,iz)];


            H = zeros(NPW, NPW); % Initialize
            PW_G = [0 0 0; Gvecs]; % Plane waves (include G=0)

            % Plane wave Hamiltonian
            for p = 1:NPW
                kp = k + PW_G(p,:);
                H(p,p) = hbar^2 * dot(kp,kp)/(2*m); % KE

                for q = 1:NPW
                    dq = PW_G(p,:) - PW_G(q,:);
                    % potential couples plane waves differing by G
                    for g = 1:size(Gvecs,1)
                        if all(abs(dq - Gvecs(g,:)) < 1e-12)
                            H(p,q) = H(p,q) + V_G(g);
                        end
                    end
                end
            end

            Evals = sort(eig(H));
            Ebands(ix,iy,iz) = Evals(1); % Only lowest band kept for FS

        end
    end
end

% Fermi Energy
EF = .5;
%% Plotting 3D Fermi Surface (Reduced Zone Scheme)
figure(1); clf;
p = patch(isosurface(kx, ky, kz, Ebands, EF));
set(p, 'FaceColor', 'red', 'EdgeColor', 'none');
xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)'); zlabel('k_z (rad/m)');
title('3D Fermi Surface (Electron in Periodic Potential)');
daspect([1 1 1]); camlight; lighting gouraud; box on;

%% Plotting 3D Fermi Surfaces (Periodic Zone Scheme)
figure(2); clf; hold on; axis equal; box on;
title('3D Fermi Surface in the Periodic Zone Scheme');

% Fermi surface in 1st BZ
F1 = isosurface(kx, ky, kz, Ebands, EF);
% To repeat FS
shifts = -1:1;
for nx = shifts
    for ny = shifts
        for nz = shifts
            % shift vertices by reciprocal lattice vectors, G0 = 2*pi/a
            shift_vec = [nx*G0, ny*G0, nz*G0];
            p = patch('Faces', F1.faces, ...
                'Vertices', F1.vertices + shift_vec, ...
                'FaceColor', 'red', ...
                'EdgeColor', 'none', ...
                'FaceAlpha', .7); %partly transparent
        end
    end
end

xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)'); zlabel('k_z (rad/m)');
view(3); camlight; lighting gouraud;


%% Effective Mass Tensor on FS
% Vertex Normals
function Vnorm = vertex_normals(FV)
V = FV.vertices;
F = FV.faces;
Nverts = size(V,1);
% face normals
fn = cross(V(F(:,2),:) - V(F(:,1),:), ...
    V(F(:,3),:) - V(F(:,1),:), 2);
fn = fn ./ vecnorm(fn,2,2); % normalize
% Initialize
Vnorm = zeros(Nverts,3);

% face normals to vertices
for i = 1:size(F,1)
    Vnorm(F(i,1),:) = Vnorm(F(i,1),:) + fn(i,:);
    Vnorm(F(i,2),:) = Vnorm(F(i,2),:) + fn(i,:);
    Vnorm(F(i,3),:) = Vnorm(F(i,3),:) + fn(i,:);
end
% Normalized
Vnorm = Vnorm ./ vecnorm(Vnorm,2,2);
end

Vnormals = vertex_normals(F1);

% Effective mass
hbar = 1;
dk = kx(2,1,1)-kx(1,1,1);
Einterp = griddedInterpolant(kx,ky,kz,Ebands,'spline');
Nverts = size(F1.vertices,1);
m_eff = zeros(Nverts,1);

for i = 1:Nverts
    k0 = F1.vertices(i,:);
    x=k0(1); y=k0(2); z=k0(3);
    % Discrete second derivatives of E
    Exx = (Einterp([x+dk,y,z]) - 2*Einterp([x,y,z]) + Einterp([x-dk,y,z])) / dk^2;
    Eyy = (Einterp([x,y+dk,z]) - 2*Einterp([x,y,z]) + Einterp([x,y-dk,z])) / dk^2;
    Ezz = (Einterp([x,y,z+dk]) - 2*Einterp([x,y,z]) + Einterp([x,y,z-dk])) / dk^2;

    Exy = ( Einterp([x+dk,y+dk,z]) - Einterp([x+dk,y-dk,z]) ...
        - Einterp([x-dk,y+dk,z]) + Einterp([x-dk,y-dk,z]) ) / (4*dk^2);

    Exz = ( Einterp([x+dk,y,z+dk]) - Einterp([x+dk,y,z-dk]) ...
        - Einterp([x-dk,y,z+dk]) + Einterp([x-dk,y,z-dk]) ) / (4*dk^2);

    Eyz = ( Einterp([x,y+dk,z+dk]) - Einterp([x,y+dk,z-dk]) ...
        - Einterp([x,y-dk,z+dk]) + Einterp([x,y-dk,z-dk]) ) / (4*dk^2);

    Eh = (1/hbar^2)*[Exx Exy Exz; Exy Eyy Eyz; Exz Eyz Ezz];

    n = Vnormals(i,:);
    m_eff(i) = 1/(n*Eh*n'); % Effective mass along normals
end

%% Plot
figure(3); clf;
p = patch(F1);
set(p,'FaceVertexCData',m_eff,'FaceColor','interp','EdgeColor','none');
colorbar; title('Effective Mass on 3D Fermi Surface');
axis equal; view(3); camlight; lighting gouraud;
xlabel('k_x'); ylabel('k_y'); zlabel('k_z');

%% Effective Mass on FS in Periodic Zone Scheme
figure(4); clf; hold on; axis equal; box on;
title('Effective Mass on 3D Fermi Surface in Periodic Zone Scheme');
% To repeat FS
shifts = -1:1;
for nx = shifts
    for ny = shifts
        for nz = shifts
            % shift vertices by reciprocal lattice vectors, G0 = 2*pi/a
            shift_vec = [nx*G0, ny*G0, nz*G0];
            p = patch('Faces', F1.faces, ...
                'Vertices', F1.vertices + shift_vec, ...
                'FaceVertexCData', m_eff, ...
                'FaceColor', 'interp', ...
                'EdgeColor', 'none', ...
                'FaceAlpha', 0.7);%partly transparent
        end
    end
end

xlabel('k_x (rad/m)'); ylabel('k_y (rad/m)'); zlabel('k_z (rad/m)');
view(3); camlight; lighting gouraud;
c = colorbar; c.Label.String = "Effective Mass (multiple of m_0)";

