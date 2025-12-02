%% (1) Incidence Matrix (-1 for current going into node, +1 out of)
A = [...
 -1  0  0  1; % Edge 1 (i_1 into node 1, out of node 4)...
  1  0  0 -1; 
  1 -1  0  0; 
  1  0 -1  0; 
  0 -1  0  1; 
  0  1 -1  0; 
  0  0  1 -1];% ...Edge 7

%% (2)
R = [0.01 20 6 10 12 3 5]; % Edge resistances
C = diag(1 ./ R); % Edge Conductances
b = [-12; 0; 0; 0; 0; 0; 0]; % Edge voltage sources, edge 1 has 12 V

% Solve using f = 0, then A^T Cb - A^T CAu = 0, so A^T C A u = A^T C b,
% or u = (A^T C A) \ (A^T C b), the voltage of each node
L = A' * C * A; % A^T C A     
rhs = A' * C * b; % A^T C b

L_reduced = L(1:3,1:3); % V4 is ground, so 4th row and column not needed
rhs_reduced = rhs(1:3);
u_reduced = L_reduced \ rhs_reduced; % Solving system 
u = [u_reduced; 0]; % Add in ground for V4
disp('Node voltages u:');
disp(u);

% Edge voltages and currents
e = b - A * u;
w = C * e; % Edge currents
disp('Edge currents w:');
disp(w);

%% (3) KCL at V2
kcl_V2 = sum(A(:,2) .* w); % Sum of currents at node 2
disp('KCL for Node 2:');
disp(kcl_V2); % Essentially 0

%% (4) KVL on Edge 3 
V_3 = e(3);            % Edge voltage
IR = w(3) * R(3);   % I*R
disp('KVL for Node 3:');
disp(V_3); disp(IR);

%% (5) 1 Amp source into V2 f
f = [0; 1; 0; -1]; % 1A into V2 from V4 (ground)

% Now solving A^T C A u = A^T C b - f
rhs2 = A' * C * b - f;
rhs2_reduced = rhs2(1:3); % V4 is ground, so 4th row and column not needed
u2_reduced = L_reduced \ rhs2_reduced; % Solving system
u2 = [u2_reduced; 0];  % Add in ground for V4
disp('1 A source into V2 Node Voltages u:');
disp(u2);

% Edge voltages and currents
e2 = b - A * u2;
w2 = C * e2; % Edge currents
disp('1 A source into V2 Edge Currents w:');
disp(w2);
% KCL at V2
kcl2_V2 = sum(A(:,2) .* w2);  % Sum of currents at node 2
disp('KCL for Node 2:');
disp(kcl2_V2); % 1 A