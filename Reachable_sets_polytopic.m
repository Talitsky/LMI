clear; clc; 
% Define Uncertainty
a11 = [0.1 0.2 0.3];
a12 = [0.1 0.2 0.3];
a21 = [0.1 0.2 0.3];
a22 = [0.1 0.2 0.3];
bw1 = [0.1 0.2 0.3];
bw2 = [0.1 0.2 0.3];
bu1 = [0.1 0.2 0.3];
bu2 = [0.1 0.2 0.3];

% Define dot_x = Ax + Bw* w + Bu* u
A = [0 1; 0.0001 0];
Bw = [1; 2];
Bu = [0; 1];


Q = [10 1; 0 10];

P = [];
mat = [];
eps = 1.e-5;
Y = sdpvar(1, 2);
 

for i = 1:3 
    delta_A = [a11(i) a12(i); a21(i) a22(i)]; 
    A_i = A + delta_A;
    B_w = Bw + [bw1(i); bw2(i)];
    B_u = Bu + [bu1(i); bu2(i)];
    mat = [mat, Q*A_i' + A_i*Q + B_u* Y + Y'* B_u' + B_w*B_w'   <= 0];
end

ops = sdpsettings('solver','mosek', 'verbose', 0);
sol = optimize(mat, [], ops);
K = value(Y)*inv(Q)
