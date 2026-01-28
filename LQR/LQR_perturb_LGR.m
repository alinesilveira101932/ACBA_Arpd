%% LQR - Estimativa perturbação - aeropendulo
% USANDO LGR
% Aline Isabel - 13/01/2026
clc;clear; close all;

% ddtheta = -b*dtheta-a*sin(theta)+c*w^2
% w = sqrt((a*sen(theta)+u)/c) -> ddtheta = -b*dtheta+u
%----------------------------------------------------
% Parametros do aeropendulo
a = 3.547623e+01; b = 8.695879e-01; c = 1.097177e-06;

Ac = [0 1; 0 -b];
Bc = [0;1];
Ec = [0;1];
C = [1 0];

% Discretizacao
T = 0.1;
[A,B] = c2dm(Ac,Bc,[],[],T,'zoh');
[~,E] = c2dm(Ac,Ec,[],[],T,'zoh');

% Matriz de parametro incerto
Adc = [0 1;0 -b*1.05];
[Ad,~] = c2dm(Adc,Bc,[],[],T,'zoh');

% LGR
num = [2.359e-5 5e-5 2.5e-5 0];
den = [1 -4 6.01 -4 1];
H  = tf(num, den, T);
rlocus(H), zgrid

% LQR
Q = 10*[1 0;0 1]; R = 0.1;
q  = 15;

[K,~,~] = dlqr(A,B,q*C'*C,1);

% Ganho do observador
Achi = [A E;zeros(1,2) 1];
Bchi = [B;0]; H = eye(2);
Hchi = [H zeros(2,1)];
eig_des = [0:2]*1e-6;
L = (place(Achi',Hchi',eig_des))';
M = inv(Achi)*L;

%% Valores de equilibrio
N = inv([A - eye(2) B;C 0])*[0;0;1];
Nx = N(1:2); Nu = N(3);
aux = inv([A - eye(2) B;C 0])*[-E;0];
Mx = aux(1:2); Mu = aux(3);

%% Simulacao
x(:,1) = [0.5;0];
chi_pri(:,1) = [x(:,1);0];

% Inicializacao
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45
kend = 100;
rbar = 0;
sys = ss(Adc,Bc,[],[]);
trec = 0; xrec = x(:,1)';   nt = 10; 

for k = 1:kend
    zpri(:,k) = Hchi*chi_pri(:,k);
    z(:,k) = x(:,k);
    chi_pos(:,k) = chi_pri(:,k) + M*(z(:,k) - zpri(:,k));
    dbar_pos(:,k) = 0;
    xbar(:,k) = Nx*rbar + Mx*dbar_pos(:,k);
    ubar(:,k) = Nu*rbar + Mu*dbar_pos(:,k);
    ulqr(:,k) = ubar(:,k) - K*(x(:,k) - xbar(:,k));
    chi_pri(:,k+1) = Achi*chi_pos(:,k) + Bchi*ulqr(:,k);
    utotal(:,k) = ulqr(:,k)+a*sin(x(1,k));
    % w(:,k) = sqrt(utotal(:,k)/c);
    ulim = 3900;
    % 
    % if utotal(:,k) > ulim
    %     utotal(:,k) = ulim;
    %     ulqr(:,k) = utotal(:,k) - a*sin(x(1,k));
    % elseif utotal(:,k)<0
    %     utotal(:,k) = 0;
    %     ulqr(:,k) = utotal(:,k) - a*sin(x(1,k));
    % end

    % % simulacao linear
    % x(:,k+1) = Ad*x(:,k)+B*u(:,k);
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % t = linspace((k-1)*T,k*T,nt);
    % U = repmat(u(:,k)',nt,1);
    % [~,tout,xout] = lsim(sys,U,t,x(:,k));
    % trec = [trec;tout(2:end)];
    % xrec = [xrec;xout(2:end,:)];
    % x(:,k+1) = xrec(end,:)';

    % atualizando planta
	xini = x(:,k);
    [t, xd] = ode45(@(t,x) simulation(x, ulqr(k),a,b,c,Ac,Bc,Ec,dbar_pos(k)),[0 T], xini, options2); 
    x(:,k+1) = xd(length(t),:)';
end

% Simulacao - tempo continuo
function xdot = simulation(x,u,a,b,c,Ac,Bc,Ec,d)
w = sqrt((a*sin(x(1))+u)/c);
xdot_1 = x(2);
xdot_2 = -a*sin(x(1))-b*x(2)+c*w^2;
xdot = [xdot_1;xdot_2];
% xdot = Ac*x+Bc*u;
end

% Plots
k = 0:kend;

% x1
figure (1)
plot(k,x(1,:),LineWidth=2)
hold on; grid on; grid minor; 
plot(k,x(2,:),LineWidth=2)
legend('Theta [rad]','dTheta [rad/s]')

% % x2
% figure(2)
% plot(k,x(2,:),LineWidth=2)
% hold on; grid on; grid minor
% legend('dTheta [rad/s]')

% d
figure(6)
plot(dbar_pos,LineWidth=2)
hold on; grid on; grid minor
legend('Perturbação estimada')

% u
figure(5)
plot(ulqr,LineWidth=2)
hold on; grid on; grid minor
legend('ulqr')


% w

figure(7)
 w= abs(sqrt(utotal/c));
plot(w,LineWidth=2)
hold on; grid on; grid minor
legend('w [rpm] = sqrt(utotal/c)')

% w

figure(8)
plot(utotal,LineWidth=2)
hold on; grid on; grid minor
legend('utotal = ulqr + a*sin(theta)')

figure(9)
plot(ulqr,LineWidth=2)
hold on; grid on; grid minor
plot(a*sin(x(1,:)),LineWidth=2)
legend('ulqr','asin(theta)')