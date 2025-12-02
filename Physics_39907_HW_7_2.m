%% Data
x = [0 1 2 3 4 5]';
y = [1 21 100 240 439 701]';

%% (1 & 2) Polyfit and Normal Form Equation
A = [x.^2, x, ones(size(x))];
u_hat = (A' * A) \ (A' * y); % Normal Form solution 
p = polyfit(x, y, 2); % Polyfit
disp('Polyfit:');
disp(p');
disp('Normal Form Equation Solution:');
disp(u_hat);

%% (3) c = 0 Fit
A2 = [x.^2, x]; 
u2 = (A2' * A2) \ (A2' * y);
a2 = u2(1);
b2 = u2(2);
disp('Fit with c = 0:');
disp(a2); disp(b2); % Fit doesn't agree with theory of a = 30, b = -10, 

%% (4) c = 0, b = 0 Fit
A3 = x.^2;
a3 = (A3' * A3) \ (A3' * y);
disp('Fit with c=0, b=0:');
disp(a3); % Close to theoretical value of a = 28

%% (5) Plotting
clf;hold on;
plot(x, y, 'ko'); % Data
y1 = a2 * x.^2 + b2 * x; % Theory 1 (c = 0)
plot(x, y1, 'b-');
y2 = a3 * x.^2; % Theory 2 (c = b = 0)
plot(x, y2, 'r--');
xlabel('x'); ylabel('y');
legend('Data', 'Theory 1: c = 0', 'Theory 2: c = b = 0');
title('Polynomial Fits to Data');

