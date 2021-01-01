%% Egyszerűsítések után:
%W_P=
%
%
%
%% feladatspecifikáció

% szabályozott szakasz
T_1 = 0.1;
T_2 = 0.01;
T_3 = 0.001;

% fázistartalék
phi_t = 27.183/180*pi;

%% időállandók beállítása

% T_2 kiejtése T_D
T_D = T_2;

% szűrőparaméter
n = 0.2;

% fázismenet
syms omega T_I

phi = -pi -2*atan(T_1*omega)-atan(T_3*omega)+atan(T_I*omega)-0.1*omega;
dphi = diff(phi,omega);

sol = vpasolve([phi_t==pi+phi dphi==0],[omega T_I],[0 2*T_1*(1+sin(phi_t))/(1-sin(phi_t))+T_3*(1+sin(phi_t))/(1-sin(phi_t))]);
T_I = double(sol.T_I)
omega_c = double(sol.omega)

% erősítés megtervezése
s = j*omega_c;
P=1/abs((T_I*s+1)/(s^2*(T_1*s+1)^2*(T_3*s+1))*exp(-0.1*s))

%% ellenőrzés
s=tf('s');
W_c = P*(T_I*s+1)*(T_D*s+1)/(s*(n*T_D*s+1));
W_p = (0.002*s+1)/(s*(0.1*s+1)^2*(0.01*s+1));
W_s = 1/(0.001*s+1)*exp(-0.1*s);
W_x = W_c*W_p*W_s;
W_cl = W_c*W_p/(1+W_x);
[Gm,Pm,Wcg,Wcp] = margin(W_x);
margin(W_x)
xlim([omega_c/20 omega_c*10])
disp('Fázistartalék:')
disp(Pm)