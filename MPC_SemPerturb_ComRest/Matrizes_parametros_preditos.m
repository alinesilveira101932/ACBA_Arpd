%% Matrizes de Parametros - controle preditivo
% Auxiliar
% Autor: Aline Isabel
% Data: 03/10/2025

%% Abarra
Abarra = A;
for i = 2:N
    Abarra = [Abarra;A^i];
end

%% Bbarra

for i = 1:N
    for j = 1:N
        if i >= j
            Bbarra{i,j} = A^(i-j)*B;
        else
            Bbarra{i,j} = zeros(nx,1);
        end
    end
end
Bbarra = cell2mat(Bbarra);


%% Matrizes de peso
Rbarra = R;
for i = 2:N
    Rbarra = blkdiag(Rbarra,R);
end

Qbarra = Q;
for i = 2:N-1
    Qbarra = blkdiag(Qbarra,Q);
end
Qbarra = blkdiag(Qbarra,P); % P é Restricao final para k = N

%% Matrizes das restricoes de U
% Sxu_x_barra * x + Sxu_u_barra * u < bxu_barra

Sxu_u_barra = Vu;
for i = 2:N  
    Sxu_u_barra = blkdiag(Sxu_u_barra,Vu);
end

Sxu_x_barra = Vx;
for i = 2:N-1  
    Sxu_x_barra = blkdiag(Sxu_x_barra,Vx);
end
Sxu_x_barra = [zeros(nx,nx*N); Sxu_x_barra zeros((N-1)*nx,nx)];

%% Matrizes das restricoes de X
% Sx_barra * x < bx_barra
Sx_barra = Sx;
for i = 2:N-1
    Sx_barra = blkdiag(Sx_barra,Sx);
end

Sx_barra = blkdiag(Sx_barra,Sf); % Sf*x < bf é Restricao final para k = N






