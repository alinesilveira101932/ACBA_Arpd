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

T = 0.1;

% Dimencoes da planta
nx = length(Ac);
nu = size(Bc,2);
ny = size(C,1);
nq = size(C,1);

Td_smc = 0.1; % Tempo de amostragem da planta
N = 25; % Horizonte de predicao

% Matrizes de Peso
q = 15;
Q = q*C'*C;
R = 1;

%% Discretizacao

[A, B] = c2dm(Ac,Bc,[],[],T, 'zoh');

% Matriz de parametro incerto
Adc = [0 1;0 -b*1.05];
[Ad,~] = c2dm(Adc,Bc,[],[],T,'zoh');

% LQR
[K,P] = dlqr(A,B,Q,R);
Af = A-B*K;



%% Valores de equilibrio
aux1 = inv([A - eye(2) B;C 0])*[0;0;1];
Nx = aux1(1:2); Nu = aux1(3);

%% Restricoes Sxu *[x;u] <= bxu ; Sxu *[xbar;ubar] <= bxu
umax = 10.7;
umin = 0;

Sxu = [1 0 0; -1 0 0; a 0 1; -2*a/pi 0 -1];
bxu = [pi/2; 0; umax; -umin];
exu = 1e-3*ones(size(Sxu,1),1);

% Particionando a restricao
Vx = Sxu(3:4,1:2);
Vu = Sxu(3:4,3);
bv = bxu(3:4);

%% Restrições: Sx*x <= bx
Sx = Sxu(1:nx,1:nx); 
bx = bxu(1:nx);

