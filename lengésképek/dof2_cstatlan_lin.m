clear all
syms m ,syms l,syms k,syms F, syms w
m=1;
l=0.1;
k=5;
F=10;

M=[4/3*m*l^2 , 1/2*m*l^2; 1/2*m*l^2, 1/3*m*l^2];
K_sum=[2*k-F*l,-k+F*l; -k, k];
eqn=det(-w^2*M+K_sum)==0;
W=cast(solve(eqn,w),'double');
w1_k=W(1)   %k-mint követő
G=(-w1_k^2*M+K_sum);
A12_k=G(1,1)/G(1,2)
w2_k=W(2)
G=(-w2_k^2*M+K_sum);
A22_k=G(1,1)/G(1,2)

M=[4/3*m*l^2 , 1/2*m*l^2; 1/2*m*l^2, 1/3*m*l^2];
K_sum=[2*k-F*l,-k; -k, k-F*l];
eqn=det(-w^2*M+K_sum)==0;
W=cast(solve(eqn,w),'double');
w1_n=W(1) %n mint nem követő
G=(-w1_n^2*M+K_sum);
A12_n=G(1,1)/G(1,2)
w2_n=W(2)
G=(-w2_n^2*M+K_sum);
A22_n=G(1,1)/G(1,2)

%ábrázolás
a=0.2;
X1_k=[l*cos(1*a)+l*cos(A12_k*a),l*cos(1*a),0,l*cos(1*a),l*cos(1*a)+l*cos(A12_k*a)]
Y1_k=[l*sin(1*a)+l*sin(A12_k*a),l*sin(1*a),0,-l*sin(1*a),-l*sin(1*a)-l*sin(A12_k*a)]

X2_k=[l*cos(1*a)+l*cos(A22_k*a),l*cos(1*a),0,l*cos(1*a),l*cos(1*a)+l*cos(A22_k*a)]
Y2_k=[l*sin(1*a)+l*sin(A22_k*a),l*sin(1*a),0,-l*sin(1*a),-l*sin(1*a)-l*sin(A22_k*a)]
figure(1);
plot(X1_k, Y1_k, X2_k, Y2_k,'LineWidth',2)
axis equal 
title('Követő erővel gerjesztett rendszer lengésképei')

X1_n=[l*cos(1*a)+l*cos(A12_n*a),l*cos(1*a),0,l*cos(1*a),l*cos(1*a)+l*cos(A12_n*a)]
Y1_n=[l*sin(1*a)+l*sin(A12_n*a),l*sin(1*a),0,-l*sin(1*a),-l*sin(1*a)-l*sin(A12_n*a)]

X2_n=[l*cos(1*a)+l*cos(A22_n*a),l*cos(1*a),0,l*cos(1*a),l*cos(1*a)+l*cos(A22_n*a)]
Y2_n=[l*sin(1*a)+l*sin(A22_n*a),l*sin(1*a),0,-l*sin(1*a),-l*sin(1*a)-l*sin(A22_n*a)]

figure(2);
plot(X1_n, Y1_n,X2_n, Y2_n,'LineWidth',2)
axis equal 
title('Fix erővel gerjesztett rendszer lengésképei')

