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

%% Experimental Parmeters
N=80;  % number of particles
D=2;   % diameter of particles
Dn=D*(1+0.1*randn(1,N)); % Poly-dispersity of 10% about D
%%Dn=D*ones(1,N);
Dmax = max(Dn); % Maximum of Dn array


K=100; % spring constant for harmonic force law
rho = 1; % Areal density
%M=3;   % mass of particles
M = rho*pi*Dn.^2/4; % Polydisperse mass
E = 2; % Initial Total Energy

Lx=9*Dmax;  % size of box
Ly=11*Dmax;

phimono = (pi*N*D^2)/(4*Lx*Ly); % Monodispersity packing fraction using D = 2
%phipoly1 = (pi*(sum(Dn.^2)))/(4*Lx*Ly) % Polydispersity phi, using sum of
%Dn's where Dn = D*ones(1,N), which reduces to monodisperse case


% For the Area Fraction:
% The summation over n reduces to N * pi * D^2, since D_n = D (mono-disperse)
% This results in the following formula for the area fraction:
% phi (area fraction) = 80 * pi * (2^2) / (4 * 18 * 22)

TT=100; % total simulation time

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

    ii=y<Dn/2;
    dw=y(ii)-Dn(ii)/2;  % Bottom wall
    Fy(ii)=Fy(ii)-K*dw;
    Ewp(nt,2)=K*sum(dw.^2)/2; %PE

    ii=x>Lx-Dn/2;
    dw=x(ii)-(Lx-Dn(ii)/2);  % Right wall
    Fx(ii)=Fx(ii)-K*dw;
    Ewp(nt,3)=K*sum(dw.^2)/2; %PE

    ii=y>Ly-Dn/2;
    dw=y(ii)-(Ly-Dn(ii)/2);  % Top wall
    Fy(ii)=Fy(ii)-K*dw;
    Ewp(nt,4)=K*sum(dw.^2)/2; %PE

    ax=Fx./M;
    ay=Fy./M-g;

    vx=vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
    vy=vy+(ay_old+ay).*dt/2;

    % Kinetic energy
    Ek(nt)= sum(M .* (vx.^2+vy.^2))/2;

    ax_old=ax;
    ay_old=ay;
end


time = dt * (0:Nt-1); % Simulation Time

%%% Plot of PE, KE, and Total Energy

clf;
plot(time, Ek, 'b-'); % Kinetic Energy
hold on;
plot(time, Ep + sum(Ewp, 2), 'g-'); % Potential Energy (Total)
plot(time, Ek + Ep + sum(Ewp, 2), 'r-'); % Total energy
xlabel('Time');
ylabel('Energy');
legend('Kinetic Energy', 'Potential Energy', 'Total Energy');
title('Energy vs Time');


%%% Plot of Total Energy

% clf;
% plot(time, Ek + Ep + sum(Ewp, 2), 'r-'); % Total energy
% xlabel('Time');
% ylabel('Total Energy');
% title('Total Energy vs Time'); % Not conserved exactly due to numerical errors



