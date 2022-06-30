A=[0 1;5 -6];
b=[4;4];
C=[2 3];

U=ctrb(A,b);
r1=rank(U)
Q=obsv(A,C);
r2=rank(Q)