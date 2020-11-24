clear; clc; 
% Initial variables
a11 = [1 2 3];
a12 = [2 5 4];
a21 = [4 5 4];
a22 = [6 4 7];

b1 = [0 1];
b2 = [2 3]; 
% Define variables
G = sdpvar(2,2 , 'full');
L = sdpvar( 1, 2);
P = [];
mat = [];
eps = 1.e-5
for i = 1:3
    for j = 1:2
        A = [a11(i) a12(i); a21(i) a22(i)];
        B = [b1(j); b2(j)];
        Pij = sdpvar(2, 2);
        P = [P, Pij];
        mat = [mat, Pij>=eps*eye(2), [Pij A*G - B*L; G'*A' - L'*B'  G + G' + Pij] <= eps*eye(4)];
    end
end

ops = sdpsettings('solver','sedumi' );
sol = optimize(mat, []);
%The output is:
K = -value(L)*inv(value(G))
