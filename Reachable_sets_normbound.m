clear; clc; 

% Define dot_x = Ax + Bw* w + Bu* u
A = [-100 1; 0.0001 -100];
Bw = [1; 0];
Bu = [1; 0];
Bp = [1; 1];
Dqu = zeros(2, 1);
Cq  = [0 1];

Q = [12 0; 0 10];

P = [];
mat = [];
eps = 1.e-5;
mu = sdpvar(1, 1);
sig  = sdpvar(1, 1); 

constraints = [mu >= 0];
constraints = [constraints, [Q*A'  + A*Q + sig*Bu*Bu' + Bw*Bw'+ mu*Bp*Bp'  (Cq*Q + mu*Dqu*Bp'- sig*Dqu*Bu')';
Cq*Q + mu*Dqu*Bp'- sig*Dqu*Bu'  -mu*eye(2)]<= 0];
ops = sdpsettings('solver', 'mosek', 'verbose', 0);
sol = optimize(constraints, [], ops);
K = -0.5*value(sig)*Bu'*Q
