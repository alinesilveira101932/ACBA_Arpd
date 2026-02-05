%% Controle Preditivo - Aline
% Codigo geral
% 03/10/2025
%-----------------------------------
clc; clearvars -except Oinf1; close all;

% Load Parametros da planta
Parametros

%% Determinar O_psi_inf - todas as refs

% Spsi = [x u xbar ubar]
Spsi = blkdiag(Sxu,Sxu);

bpsi = [bxu ;bxu - exu];

Gamma = [
    eye(nx)        zeros(nx,nq)        zeros(nx,nd);          % x
    -K        (K*Nx+Nu)        (K*Mx+Mu);                       % u
    zeros(nx)        Nx Mx;           % xbar
    zeros(nu,nx)        Nu                 Mu                  % ubar
];
max_iter = 50;

Apsi = [Af B*(K*Nx+Nu) B*(K*Mx+Mu)+E; zeros(nq,nx) eye(nq) zeros(nq,nd); zeros(nd,nx) zeros(nd,nq) eye(nd)];

[So,bo,Si,bi] = Oinf_MAS(Apsi,Gamma,Spsi,bpsi,max_iter);

Ox = So(:,1:nx);
Or = So(:,nx+1);
Od = So(:,nx+2:end);

O_psi_inf = Polyhedron('H',[So bo]); % Declaração de Oinf como objeto do MPT Toolbox

Sf = Ox;

% Teste rbar = 1*pi/180 e dbar = 0;
dbar_teste = 0; rbar_teste = 1*pi/180;
% xbar_teste = Nx*rbar_teste + Mx*dbar_teste;
% ubar_teste = Nu*rbar_teste + Mu*dbar_teste;
bf_teste = bo-Or*rbar_teste-Od*dbar_teste;
O_inf = Polyhedron('H',[Sf bf_teste]); % Declaração de Oinf como objeto do MPT Toolbox
figure (120); plot(O_inf); xlabel('x1');ylabel('x2'); title('rbar = 1 [deg]; dbar = 0')

%% Dual Mode Predictive Control

% Carregando matrizes de parametros do controle preditivo
Matrizes_parametros_preditos

% Configurando solvers
options2 = odeset('Reltol',1e-6,'AbsTol',1e-6); %ode45
options_qudprog = optimoptions('quadprog','Display','off');

% Inicializacao
kend = 100;
x0 = [1;0];
x(:,1) = x0; 
chi_pri(:,1) = [x(:,1);0];
rbar = 1*pi/180;

for k = 1:kend

  % Matrizes finais do quadprog 
    % J = f*U'+0.5*U'*H*U
    % s.a A * U < B

    zpri(:,k) = Hchi*chi_pri(:,k);
    z(:,k) = x(:,k);
    chi_pos(:,k) = chi_pri(:,k) + M*(z(:,k) - zpri(:,k));
    dbar_pos(:,k) = chi_pos(3,k);
    xbar(:,k) = Nx*rbar + Mx*dbar_pos(:,k);
    ubar(:,k) = Nu*rbar + Mu*dbar_pos(:,k);

    bf = bo-Or*rbar-Od*dbar_pos(:,k);

    bxu_barra = [bv-Vx*x(:,k); repmat(bv,N-1,1)];
    bx_barra = [repmat(bx,N-1,1);bf];

    Aqp_xu = Sxu_u_barra + Sxu_x_barra*Bbarra;
    bqp_xu = bxu_barra- Sxu_x_barra*(Abarra*x(:,k)+En*repmat(dbar_pos(k),N,1));

    Aqp_x = Sx_barra*Bbarra;
    bqp_x = bx_barra- Sx_barra*(Abarra*x(:,k)+En*repmat(dbar_pos(k),N,1));

    Aqp = [Aqp_xu; Aqp_x];
    bqp = [bqp_xu;bqp_x];

    Hqp = Bbarra'*Qbarra*Bbarra + Rbarra;

    % Simetrização
    Hqp = (Hqp + Hqp')/2;

    fqp = (Bbarra'*Qbarra*(Abarra*x(:,k)+En*repmat(dbar_pos(k),N,1) - repmat(xbar(:,k),N,1)) - Rbarra*repmat(ubar(:,k),N,1));

    dummy = quadprog(Hqp, fqp,Aqp,bqp, [], [], [], [], [], options_qudprog);

    umpc(k) = dummy(1);
    utotal(k) = umpc(k) + a*sin(x(1,k));
    
    chi_pri(:,k+1) = Achi*chi_pos(:,k) + Bchi*umpc(:,k);

	xini = x(:,k);
    [t, xd] = ode45(@(t,x) simulation(x, umpc(k), a,b,c,Ac,Bc),[0 Td_smc], xini, options2); 
    x(:,k+1) = xd(length(t),:)';

end

% Simulacao - tempo continuo
function xdot = simulation(x,u,a,b,c,Ac,Bc)
w = sqrt((a*sin(x(1))+u)/c);
xdot_1 = x(2);
xdot_2 = -a*sin(x(1))-b*x(2)+c*w^2;
xdot = [xdot_1;xdot_2];
% xdot = 1.05*Ac*x+Bc*u;
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
plot(x(1,:)*180/pi,'b','LineWidth',2)
hold on;
plot(x(2,:)*180/pi,'r','LineWidth',2)
xlabel('k')
legend('x1 [deg]','x2 [deg/s]')
xlim([0 100])
% Plot do limite
hold on;
plot([0 100],[rbar*180/pi rbar*180/pi],'--k','linewidth',1,'handlevisibility','off')
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
ylim([-0.2 umax+3])
plot([0 100],[umax umax],'--r','linewidth',1,'handlevisibility','off')
text(0.1,umax+0.7,'Limitante de controle','Color','r','FontWeight','bold')
plot([0 100],[0 0],'--r','linewidth',1,'handlevisibility','off')

figure(9)
plot(umpc,LineWidth=2)
hold on; grid on; grid minor
plot(a*sin(x(1,:)),LineWidth=2)
legend('umpc','asin(theta)')
xlim([0 100])
xlabel('k')

figure(10)
plot(dbar_pos,LineWidth=2)
legend('dbar');grid on;grid minor
xlim([0 100])
xlabel('k')