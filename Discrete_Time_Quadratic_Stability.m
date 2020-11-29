clear; clc; 
 
a11 = [0.1 0.2 0.3];
a12 = [0.1 0.2 0.3];
a21 = [0.1 0.2 0.3];
a22 = [0.1 0.2 0.3];
A = [0 1; 0.0001 0];

P = [];
mat = [];
eps = 1.e-5;
P = sdpvar(2, 2);
mat = [mat, P >= eps*eye(2)];

for i = 1:3 
    delta_A = [a11(i) a12(i); a21(i) a22(i)]; 
 
    mat = [mat,   ((A + delta_A)') * P * (A + delta_A) - P <= 0];
 
end

ops = sdpsettings('solver','mosek', 'verbose', 0);
sol = optimize(mat, [], ops);
P = value(P) 
