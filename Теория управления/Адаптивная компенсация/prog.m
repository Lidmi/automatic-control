A=[2 8;3 6];
b=[2;3];
U=ctrb(A,b);
rank(U);
Ad=[0 1;-256 -32];
H=[1 0];
U2=obsv(Ad,H);
rank(U2);

M=sylv(A,-Ad,b*H)
K=H*inv(M)

Af0=[0 1; -0.01 -sqrt(0.02)];
bf0=[0;1];
syms n11 n12 n21 n22;
n=[n11 n12;n21 n22];
eqn=n*b==bf0;
sol=solve([eqn],[n11, n12, n21, n22]);
N=[double(sol.n11) double(sol.n12); double(sol.n21) double(sol.n22)];

Am=A-b*K;
Q=eye(2);
P=lyap(Am',Q);

