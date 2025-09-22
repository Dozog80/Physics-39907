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
K=100; % spring constant for harmonic force law
M=3;   % mass of particles

Lx=9*D;  % size of box
Ly=11*D;

% For the Area Fraction:
% The summation over n reduces to N * pi * D^2, since D_n = D (mono-disperse)
% This results in the following formula for the area fraction:
% phi (area fraction) = 80 * pi * (2^2) / (4 * 18 * 22)

Arf = N * pi * D^2 / (4 * Lx * Ly)

Nx = Lx/D; % Number of particles along the x-direction
Ny = Ly/D; % Number of particles along the y-direction
Nmax = Nx * Ny; % Maximum number of particles in rectangular container
Arfmax = Nmax * pi * D^2 / (4 * Lx * Ly) % Max Area Fraction = 0.7854

TT=100; % total simulation time

%% Physical Parameters
g=0;

%% Display Parameters
plotit=true;  % plot ?
Nplotskip=50;  % number of timesteps to skip before plotting

%% Simulation Parmeters
dt=1e-2;
Nt=fix(TT/dt); % Number of time steps

%% Initial Conditions
[x y]=ndgrid(D/2:D:Lx-D/2,D/2:D:Ly-D/2);
ii=randperm(numel(x),N);
x=x(ii(1:N));
y=y(ii(1:N));

vx=randn(1,N)/3;
vy=randn(1,N)/3;

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
    h(np)=rectangle('Position',[x(np)-.5*D y(np)-.5*D D D],'Curvature',[1 1],'edgecolor','b');
  end
  axis('equal');
  axis([0 Lx 0 Ly]);
end

%% Main Loop

for nt=1:Nt
  
  % plot particles
  if(plotit && rem(nt-1,Nplotskip)==0)
    for np=1:N
      set(h(np),'Position',[x(np)-.5*D y(np)-.5*D D D]);
    end
    drawnow;
    % If we add the area fraction formula (Arf = N * pi * D^2 / (4 * Lx *
    % Ly)) here we see that the area fraction does not change during the
    % simulation. Since the particle size and the area of the container
    % does not change during the simulation, and particles are not leaving 
    % or entering the containter, there will be no change in the area 
    % fraction. The particles also don't overlap with each other or extend
    % beyond the walls of the container, so during the simulation the area
    % the particles cover doesn't change. 
  end
  
  x=x+vx*dt+ax_old.*dt.^2/2;  % first step in Verlet integration
  y=y+vy*dt+ay_old.*dt.^2/2;
  
  % Interaction detector and Force Law
  Fx=zeros(1,N);
  Fy=zeros(1,N);
  
  for nn=1:N
    for mm=nn+1:N
      dy=y(mm)-y(nn);
      if(abs(dy)<=D)
        dx=x(mm)-x(nn);
        dnm=dx.^2+dy.^2;
        if(dnm<D^2)
          dnm=sqrt(dnm);
          F=-K*(D/dnm-1);
          Ep(nt)=Ep(nt)+(D-dnm).^2;  % particle-particle PE
          Fx(nn)=Fx(nn)+F.*dx;  % particle-particle Force Law
          Fx(mm)=Fx(mm)-F.*dx;
          Fy(nn)=Fy(nn)+F.*dy;  % particle-particle Force Law
          Fy(mm)=Fy(mm)-F.*dy;
          
        end
      end
    end
  end
  
  Ep(nt)=K/2*Ep(nt);
  
  ii=x<D/2;
  dw=x(ii)-D/2; % Left wall
  Fx(ii)=Fx(ii)-K*dw;
  Ewp(nt,1)=K*sum(dw.^2)/2; %PE
  
  ii=y<D/2;
  dw=y(ii)-D/2;  % Bottom wall
  Fy(ii)=Fy(ii)-K*dw;
  Ewp(nt,2)=K*sum(dw.^2)/2; %PE
  
  ii=x>Lx-D/2;
  dw=x(ii)-(Lx-D/2);  % Right wall
  Fx(ii)=Fx(ii)-K*dw;
  Ewp(nt,3)=K*sum(dw.^2)/2; %PE
  
  ii=y>Ly-D/2;
  dw=y(ii)-(Ly-D/2);  % Top wall
  Fy(ii)=Fy(ii)-K*dw;
  Ewp(nt,4)=K*sum(dw.^2)/2; %PE
  
  ax=Fx./M;
  ay=Fy./M-g;
  
  vx=vx+(ax_old+ax).*dt/2;  % second step in Verlet integration
  vy=vy+(ay_old+ay).*dt/2;
  
  % Kinetic energy
  Ek(nt)=M*sum((vx.^2+vy.^2))/2;
  
  ax_old=ax;
  ay_old=ay;
end


%%% Question 1
% 1)
% The summation over n reduces to N * pi * D^2, since D_n = D (mono-disperse)
% This results in the following formula for the area fraction:
% phi (area fraction) = 80 * pi * (2^2) / (4 * 18 * 22)
% Entering the right side of the equal sign into the command window results
% in an area fraction of 0.6347, which was intended to be shown.
% Code was added to calculate the area fraction for the mono-disperse case
% ^ (See Arf) ^

% 2)
% The area fraction does not change during the simulation. In the
% definition of the area fraction, the relavant experimental parameters are
% the number of particles (N), the diameter (D), and the side lengths of
% the container (Lx and Ly). None of these parameters change during the
% simulation, the numbers of particles, particle size, and rectangular
% container size are all constant so the area fraction must be constant. We
% can confirm this by adding the area fraction formula (Arf = N * pi * D^2
% / (4 * Lx * Ly)) here and seeing that the value of the area fraction does
% not change.

% 3)
% The largest area fraction we can achieve using the square packing initial
% condition can be calculated using the area fraction formula for the 
% maximum number of particles that can fit in the container, or, 
% Nmax * pi * (D^2) / (4 * Lx * Ly). Assuming the experimental parameters
% used in part 1 of this question (Lx = 9*D, Ly = 11*D, and D = 2), the
% formula reduces to 99 * pi * (2^2) / (4 * 18 * 22), where the 99 maximum
% particles is calculated by multiplying the number of particles that can
% fit along the x-direction by the number of particles that can fit along
% the y-direction (9 by 11). Code for this maximum area fraction is
% included in the experimental parameters section (^ See Arfmax ^). 

% 4) 
% The maximum realizable area fraction does depend on the shape and size of
% the container beacuse if the dimensions of the container aren't an
% integer multiple of particle diameter (if we had say a rectangle 9.5D by
% 11.5D) then there will always be some wasted space that the particles
% cannot be packed into which increases the denominator of the area
% fraction formula and results in a smaller maximum area fraction. 

% 5) 
% The theoretical maximum packing fraction is the area of particles inside
% a hexagonal primitive lattice cell divided by the area of the hexagon.
% For a hexagon of side length D, equal to the diameter of the particles,
% the packing fraction becomes (since 3 total particles fit in every 
% hexagon, one in the center, and 1/3 of a particle on every side) 
% Arfmax = (3 * pi * D^2 / 4) / (3 * sqrt(3) * D^2 / 2), where the
% denominator is the area of the hexagon. This reduces to pi / 2 * sqrt(3),
% which is equal to 0.90689968211. Achieving this theoretical maximum
% packing fraction of hexagonal packing in a finite simulation would likely
% not be possible with the large number of particles and interactions
% needed being computationally expensive, so only approximations of this
% theoretical value would be reached.
