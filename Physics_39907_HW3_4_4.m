%% Molecular Dynamics Simulator
% <md.m> Mark D. Shattuck 7/22/2010

% revision history:
% 7/22/2010 Mark D. Shattuck <mds> md.m
%           MD Demo for HandsOn 2010
%           000 mds Initial conditions and visualization
% 7/24/2010 001 mds Add Euler
% 7/24/2010 002 mds Add Interaction detection and Force Law
% 7/24/2010 003 mds Add KE, PE, Nplotskip, and plotit
% 1/21/2011 010 mds Add Velocity Verlet

phi_list = 0.02 + (0:19) / 19 * (0.7 - 0.02);  % 20 values from 0.02 to 0.7
Nphi = length(phi_list);
Zs = zeros(1, Nphi);  % preallocate results

for nphi=1:Nphi % loop over all area fractions
    phi=phi_list(nphi); % set phi to np−th value
    disp(nphi); % user feedback

    %% Experimental Parmeters
    N=80;  % number of particles
    D=2;   % diameter of particles
    pd=0; % particle diameter poly−dispersity
    Dn=D*(1+pd*randn(1,N)); % Poly-dispersity of 10% about D
    %%Dn=D*ones(1,N);
    Dmax = max(Dn); % Maximum of Dn array


    K=1000; % spring constant for harmonic force law
    rho = 1; % Areal density
    %%M=3;   % mass of particles
    M = rho*pi*Dn.^2/4; % Polydisperse mass

    % Lx=3*9*D;  % size of box
    % Ly=3*11*D;
    G = 9 / 11;  % Fixed aspect ratio Lx / Ly
    Lx = sqrt(pi * G * sum(Dn.^2) / phi) / 2;
    Ly = Lx / G;

    TT=1000; % total simulation time

    E=10; % Total Energy

    % phimono = (pi*N*D^2)/(4*Lx*Ly); % Monodispersity packing fraction using D = 2
    %phipoly1 = (pi*(sum(Dn.^2)))/(4*Lx*Ly) % Polydispersity phi, using sum of
    %Dn's where Dn = D*ones(1,N), which reduces to monodisperse case


    %% Physical Parameters
    g=0;

    %% Display Parameters
    plotit=false;  % plot ?
    Nplotskip=50;  % number of timesteps to skip before plotting

    %% Simulation Parmeters
    dt=1e-2;
    Nt=fix(TT/dt); % Number of time steps

    %% Initial Conditions
    [x y]=ndgrid(Dmax/2:Dmax:Lx-Dmax/2,Dmax/2:Dmax:Ly-Dmax/2);
    ii=randperm(numel(x),N);
    x=x(ii(1:N));
    y=y(ii(1:N));

    vx=randn(1,N);
    vy=randn(1,N);

    K0 = E; % Initially all Energy is Kinetic, since there is no overlaps
    % of the particles with the walls or each other initially, so there is no
    % potential energy initially.

    v0 = sqrt(2 * K0/sum(M)); % Initial velocity based on classical equation,
    % (K = Mv^2/2)

    % Initial Kinetic Energy of particles, not corrected to be equal to E
    % since the initial particle velocities are random
    K_initial = 0.5 * sum(M .* (vx.^2 + vy.^2));

    % Kinetic correction (kc) so that the initial particle kinetic energy
    % will be E
    kc = sqrt(K0 / K_initial);

    vx = kc * vx;
    vy = kc * vy;

    ax_old=0*x;    % inital condition Fnm=0;
    ay_old=0*y-g;

    %% Output
    Ek=zeros(Nt,1);    % Kinetic Energy
    Ep=zeros(Nt,1);    % particle-particle potential
    Ewp=zeros(Nt,4);   % wall-particle potential (1234)=>(LBRT)
    Fws = zeros(Nt, 4);  % Initialize array to hold force measurements on walls

    %% Setup Plotting
    if(plotit)
        clf;
        h=zeros(1,N);
        for np=1:N
            h(np)=rectangle('Position',[x(np)-.5*Dn y(np)-.5*Dn Dn Dn],'Curvature',[1 1],'edgecolor','b');
        end
        axis('equal');
        axis([0 Lx 0 Ly]);
    end

    %% Main Loop

    for nt=1:Nt

        % plot particles
        if(plotit && rem(nt-1,Nplotskip)==0)
            for np=1:N
                set(h(np),'Position',[x(np)-.5*Dn y(np)-.5*Dn Dn Dn]);
            end
            drawnow;
        end

        x=x+vx*dt+ax_old.*dt.^2/2;  % first step in Verlet integration
        y=y+vy*dt+ay_old.*dt.^2/2;

        % Interaction detector and Force Law
        Fx=zeros(1,N);
        Fy=zeros(1,N);


        for nn=1:N
            for mm=nn+1:N
                dy=y(mm)-y(nn);
                Dnm = ((Dn(nn) + Dn(mm))/2);
                if(abs(dy)<=(Dnm))
                    dx=x(mm)-x(nn);
                    dnm=dx.^2+dy.^2;
                    if(dnm<(Dnm)^2)
                        dnm=sqrt(dnm);
                        F=-K*(Dnm/dnm-1);
                        Ep(nt)=Ep(nt)+((Dnm)-dnm).^2;  % particle-particle PE
                        Fx(nn)=Fx(nn)+F.*dx;  % particle-particle Force Law
                        Fx(mm)=Fx(mm)-F.*dx;
                        Fy(nn)=Fy(nn)+F.*dy;  % particle-particle Force Law
                        Fy(mm)=Fy(mm)-F.*dy;
                    end
                end
            end
        end
        Ep(nt)=K/2*Ep(nt);

        ii=x<Dn/2;
        dw=x(ii)-Dn(ii)/2; % Left wall
        Fx(ii)=Fx(ii)-K*dw;
        Ewp(nt,1)=K*sum(dw.^2)/2; %PE
        Fws(nt,1) = -K * sum(dw); % left wall force from F = -dU/dr on
        % PE equation (1/2Kx^2)


        ii=y<Dn/2;
        dw=y(ii)-Dn(ii)/2;  % Bottom wall
        Fy(ii)=Fy(ii)-K*dw;
        Ewp(nt,2)=K*sum(dw.^2)/2; %PE
        Fws(nt,2) = -K * sum(dw); % bottom wall force


        ii=x>Lx-Dn/2;
        dw=x(ii)-(Lx-Dn(ii)/2);  % Right wall
        Fx(ii)=Fx(ii)-K*dw;
        Ewp(nt,3)=K*sum(dw.^2)/2; %PE
        Fws(nt,3) = -K * sum(dw); % right wall force


        ii=y>Ly-Dn/2;
        dw=y(ii)-(Ly-Dn(ii)/2);  % Top wall
        Fy(ii)=Fy(ii)-K*dw;
        Ewp(nt,4)=K*sum(dw.^2)/2; %PE
        Fws(nt,4) = -K * sum(dw); % top wall force



        ax=Fx./M;
        ay=Fy./M-g;

        vx=vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
        vy=vy+(ay_old+ay).*dt/2;

        % Kinetic energy
        Ek(nt)= sum(M .* (vx.^2+vy.^2))/2;

        ax_old=ax;
        ay_old=ay;
    end

 

    % Average Force skipping initial 10 percent of simulation
    Fl = mean(Fws(fix(0.1 * Nt):end,1));
    Fb = mean(Fws(fix(0.1 * Nt):end,2));
    Fr = mean(Fws(fix(0.1 * Nt):end,3));
    Ft = mean(Fws(fix(0.1 * Nt):end,4));

    Pl = Fl / Ly;
    Pb = Fb / Lx;
    Pr = Fr / Ly;
    Pt = Ft / Lx;

    % Average pressure across all four walls, we take mean(abs(...)) for the average pressure since the pressures should not take
    % into account the signs of the forces, just the magnitudes of the forces,

    P = mean(abs([Pl, Pr, Pt, Pb]));

    K_avg = mean(Ek(fix(0.1*end):end));  % average kinetic energy
    T = K_avg / N; % Temperature per particle, using classical definition

    A = Lx * Ly;
    Zs(nphi) = (P * A) / (N * T);
    % plot results as they are calculated
    figure(2);
    plot(phi_list,Zs,'.');
    drawnow;
    figure(1);
end

Z_mean = mean(Zs); % Running the code gave 1.1950 = the average value of Z
Z_std = std(Zs); % And 0.0144 for the standard deviation of Z



%%% (1)
% (a) Physically, the particles are always exerting an outward force on the
% walls, so the tensor convention for pressure being positive whenever the
% particles are pushing the walls outward means the pressure will always be
% positive.
% (b) The pressures differ for each wall because the system has a small
% number of particles and a finite amount of simulation time, so there will
% always be a natural and random fluctuation over time. Also the box isn't
% square, and initial velocities and positions are random, so the system
% does not have symmetry. The larger vertical size of the box is reflected
% in the slightly larger pressure on the top and bottom walls.
% (c) The numbers are not the same everytime the code is run. The random
% initial conditions and the natural fluctuations from these random
% conditions and small simulation time and particle number will result in
% slight differences in the pressures measured.

%%% (2)
% Running the code gave 1.1950 = the average value of Z
% and 0.0144 for the standard deviation of Z

%%% (5) 
% For the correction of z = 1/(1-2phi): The excluded area will be pi * D^2 
% = pi, but for each particle we half this area to get pi/2. So for N 
% particles, the effective area will be A - N * pi/2. The area packing 
% fraction will be N * pi/4 / A. Dividing the effective area by A gives 
% 1 - N * pi/2 / A = 1 - 2phi. Then z = A / Aeff = 1 / 1 - 2phi.
% Run the following code to plot the measured and theoretical values
%%% of z
phi_ideal = linspace(0, 0.75, 100);
Z_ideal = ones(size(phi_ideal));
Z_vdw = 1 ./ (1 - 2 * phi_ideal);
Z_vdw(phi_ideal >= 0.5) = NaN;  % avoid phi = 0.5
Z_linear = 1 + 2 * phi_ideal;
Z_adhoc = 1 + (2 * phi_ideal) ./ ((1 - phi_ideal).^2);
Z_cs = 1 + (2 * phi_ideal - (7/8) * phi_ideal.^2) ./ ((1 - phi_ideal).^2);

% Plotting z vs phi
plot(phi_list, Zs, 'ko', 'DisplayName', 'Measured');
hold on;
plot(phi_ideal, Z_ideal, 'r--', 'DisplayName', 'Ideal Gas: Z = 1');
plot(phi_ideal, Z_vdw, 'g--', 'DisplayName', 'Van der Waals');
plot(phi_ideal, Z_linear, 'b--', 'DisplayName', 'Linear VDW');
plot(phi_ideal, Z_adhoc, 'c--', 'DisplayName', 'Ad-hoc VDW');
plot(phi_ideal, Z_cs, 'm--', 'DisplayName', 'Carnahan and Starling');
xlabel('\phi (Area Packing Fraction)');
ylabel('Z = PA / NT');
title('Z vs \phi (Area Packing Fraction)');
legend('Location', 'NorthWest');
axis([0 0.75 0 15]);


