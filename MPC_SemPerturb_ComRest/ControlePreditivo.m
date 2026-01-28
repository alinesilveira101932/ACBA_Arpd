%% Controle Preditivo - Aline
% Codigo geral
% 03/10/2025
%-----------------------------------
clc; clear; close all;

% Load Parametros da planta
Parametros

%% Determinar O_psi_inf - todas as refs

Gamma = [eye(nx) zeros(nx,nq); -K (K*Nx+Nu); zeros(nx) Nx; zeros(nu,nx) Nu];

Spsi = blkdiag(Sxu,Sxu);
bpsi = [bxu ;bxu - exu];

Apsi_f = [Af B*(K*Nx+Nu); zeros(nq,nx) eye(nq)];
max_iter = 80;

[So,bo,Si,bi] = Oinf_MAS(Apsi_f,Gamma,Spsi,bpsi,max_iter);

O_psi_inf = Polyhedron('H',[So bo]); % Declaração de Oinf como objeto do MPT Toolbox

% figure (1)
% plot(O_psi_inf)

%% Determinar Oinf a partir de O_psi_inf
Ox = So(:,1:nx);
Or = So(:,nx+1:end);
Sf = Ox;

%% Dual Mode Predictive Control

% Carregando matrizes de parametros do controle preditivo
Matrizes_parametros_preditos

% Configurando solvers
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45
options_qudprog = optimoptions('quadprog','Display','off');

% Inicializacao
kend = 100;
x0 = [0;0];

% Loop
x(:,1) = x0;
rbar = 0.3;
sys = ss(Adc,Bc,[],[]);
trec = 0; xrec = x(:,1)'; nt = 10; 
for k = 1:kend

  % Matrizes finais do quadprog 
    % J = f*U'+0.5*U'*H*U
    % s.a A * U < B
    
    xbar = Nx*rbar;
    ubar = Nu*rbar;

    bf = bo - Or*rbar;

    bxu_barra = [bv-Vx*x(:,k); repmat(bv,N-1,1)];
    bx_barra = [repmat(bx,N-1,1);bf];

    Aqp_xu = Sxu_u_barra + Sxu_x_barra*Bbarra;
    bqp_xu = bxu_barra- Sxu_x_barra*Abarra*x(:,k);

    % Aqp_x = Sx_barra*Bbarra;
    % bqp_x = bx_barra- Sx_barra*Abarra*x(:,k);
    Aqp_x = Sf*Bbarra(end-nx+1:end,:);
    bqp_x = bf- Sf*Abarra(end-nx+1:end,:)*x(:,k);
    Aqp = [Aqp_xu; Aqp_x];
    bqp = [bqp_xu;bqp_x];

    Hqp = Bbarra'*Qbarra*Bbarra + Rbarra;

    % Simetrização
    Hqp = (Hqp + Hqp')/2;

    fqp = (Bbarra'*Qbarra*(Abarra*x(:,k) - repmat(xbar,N,1)) - Rbarra*repmat(ubar,N,1));

    dummy = quadprog(Hqp, fqp,Aqp,bqp, [], [], [], [], [], options_qudprog);

    umpc(k) = dummy(1);
    utotal(k) = umpc(k) + a*sin(x(1,k));
    w(k) = sqrt(utotal(k)/c);
% % simulacao linear
    x(:,k+1) = A*x(:,k)+B*umpc(:,k);
    t = linspace((k-1)*T,k*T,nt);
    U = repmat(umpc(:,k)',nt,1);
    [~,tout,xout] = lsim(sys,U,t,x(:,k));
    trec = [trec;tout(2:end)];
    xrec = [xrec;xout(2:end,:)];
    x(:,k+1) = xrec(end,:)';
    % Simulacao continua
	% xini = x(:,k);
    % [t, xd] = ode45(@(t,x) simulation(x, ulqr(k), a,b,c,Ac,Bc),[0 Td_smc], xini, options2); 
    % x(:,k+1) = xd(length(t),:)';

end

% Simulacao - tempo continuo
function xdot = simulation(x,u,a,b,c,Ac,Bc)
% w = sqrt((a*sin(x(1))+u)/c);
% xdot_1 = x(2);
% xdot_2 = -a*sin(x(1))-b*x(2)+c*1.05*w^2;
% xdot = [xdot_1;xdot_2];
xdot = Ac*x+Bc*u;
end

% Plots
figure (1)
hold on
grid on; 
plot(umpc,'r','LineWidth',2)
hold on;
xlabel('k')
legend('umpc')
xlim([0 100])



figure (2)
hold on
grid on; 
plot(x(1,:),'b','LineWidth',2)
hold on;
plot(x(2,:),'r','LineWidth',2)
xlabel('k')
legend('x1 [rad]','x2 [rad/s]')
xlim([0 100])
% Plot do limite
% hold on;
% plot([0 100],[0.15 0.15],'--k','linewidth',1,'handlevisibility','off')
% text(0.1,0.18,'Limitante de x2','Color','r','FontWeight','bold')

figure (4)
w = sqrt(utotal/c);
hold on
grid on; 
plot(w,'r','LineWidth',2)
hold on;
xlabel('k')
legend('w [rpm]')
xlim([0 100])

figure(8)
plot(utotal,LineWidth=2)
hold on; grid on; grid minor; xlim([0 100]);
legend('utotal = umpc + a*sin(theta)')
xlabel('k')
% Plot do limite
hold on;
ylim([0 12])
plot([0 100],[10.7 10.7],'--k','linewidth',1,'handlevisibility','off')
text(0.1,11.3,'Limitante de controle','Color','r','FontWeight','bold')


figure(9)
plot(umpc,LineWidth=2)
hold on; grid on; grid minor
plot(a*sin(x(1,:)),LineWidth=2)
legend('umpc','asin(theta)')
xlim([0 100])
xlabel('k')