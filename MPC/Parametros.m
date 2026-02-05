%% Parametros - controle preditivo - Mudanca ref + estimativa perturbacao
% Auxiliar
% AEROPENDULO
% Autor: Aline Isabel
% Data: 16/01/2026

%% Parâmetros planta

% ddtheta = -b*dtheta-a*sin(theta)+c*w^2
% w = sqrt((a*sen(theta)+u)/c) -> ddtheta = -b*dtheta+u
%----------------------------------------------------
% Parametros do aeropendulo
a = 3.547623e+01; b = 8.695879e-01; c = 1.097177e-06;

Ac = [0 1; 0 -b];
Bc = [0;1];
C = [1 0];
Ec = [0;1];
T = 0.1;

% Dimencoes da planta
nx = length(Ac);
nu = size(Bc,2);
ny = size(C,1);
nq = size(C,1);
nd = size(Ec,2);

Td_smc = 0.1; % Tempo de amostragem da planta
N = 60; % Horizonte de predicao

% Matrizes de Peso
q = 120;
Q = q*C'*C;
R = 1;

%% Discretizacao

[A, B] = c2dm(Ac,Bc,[],[],T, 'zoh');
[~, E] = c2dm(Ac,Ec,[],[],T, 'zoh');

% Matriz de parametro incerto
Adc = [0 1;0 -b*1.05];
[Ad,~] = c2dm(Adc,Bc,[],[],T,'zoh');

% LQR
[K,P] = dlqr(A,B,Q,R);
Af = A-B*K;


% Ganho do observador
Achi = [A E;zeros(1,2) 1];
Bchi = [B;0]; H = eye(2);
Hchi = [H zeros(2,1)];
eig_des = [0:2]*1e-6;
L = (place(Achi',Hchi',eig_des))';
M = inv(Achi)*L;


%% Ponto de equilibrio
aux = inv([A - eye(2) B;C 0])*[0;0;1];
Nx = aux(1:2);
Nu = aux(3);
aux2 = inv([A-eye(nx) B; C zeros(nq,nu)])*[-E; zeros(nq,nd)];
Mx = aux2(1:2);
Mu = aux2(3);

%% Restricoes Sxu *[x;u] <= bxu ; Sxu *[xbar;ubar] <= bxu
umax = 30;
umin = 0;

Sxu = [1 0 0; -1 0 0; a 0 1; -2*a/pi 0 -1];
bxu = [pi/2; 0; umax; -umin];
exu = 1e-3*ones(size(Sxu,1),1);

% Particionando a restricao
Vx = Sxu(3:4,1:2);
Vu = Sxu(3:4,3);
bv = bxu(3:4);
ev = 1e-5*ones(size(bv,1),1);

%% Restrições: Sx*x <= bx
Sx = Sxu(1:nx,1:nx); 
bx = bxu(1:nx);
ex = 1e-5*ones(size(bx,1),1);
