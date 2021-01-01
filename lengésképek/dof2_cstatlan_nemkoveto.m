clear all
syms m ,syms l,syms k,syms F, syms w
m=1;
l=1;
k=5;
F=1;

M=[4/3*m*l^2 , 1/2*m*l^2; 1/2*m*l^2, 1/3*m*l^2];
K_sum=[2*k-F*l,-k; -k, k-F*l];
eqn=det(-w^2*M+K_sum)==0;
W=cast(solve(eqn,w),'double');w1_n=W(1)
G=(-w1_n^2*M+K_sum);
A12_n=G(1,1)/G(1,2)
w2_n=W(2)
G=(-w2_n^2*M+K_sum);
A22_n=G(1,1)/G(1,2)