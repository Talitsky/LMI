clear; clc; 

% Define dot_x = Ax + Bw* w + Bu* u
A = [-100 1; 0.0001 -100];
Bw = [1; 0];
Bu = [1; 0];
Bp = [1 0; 0 1];
Dqu = zeros(2, 1);
Cq  = [0 1];

Q = [12 0; 0 10];

P = [];
mat = [];
eps = 1.e-5;
M = sdpvar(2, 2);
Y  = sdpvar(1, 2); 

constraints = [M >= 0];
constraints = [constraints, [Q*A'  + A*Q + Bu*Y + Y'*Bu' + Bw*Bw'+ Bp*M*Bp'  (Cq*Q + Dqu*Y)';
Cq*Q + Dqu*Y  -M]<= 0];
ops = sdpsettings('solver', 'mosek', 'verbose', 0);
sol = optimize(constraints, [], ops);
K = value(Y)*inv(Q)
