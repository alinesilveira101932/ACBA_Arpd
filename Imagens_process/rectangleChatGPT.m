function rectangleChatGPT(y_inf, y_sup, x_centro, largura, y_linha, cor)

% Cálculo dos limites do retângulo
x_esq = x_centro - largura/2;
x_dir = x_centro + largura/2;
altura = y_sup - y_inf;

% Criar figura
% figure;
% hold on;
% axis equal;

% Desenhar retângulo
rectangle('Position', [x_esq, y_inf, largura, altura], ...
          'EdgeColor', cor, 'LineWidth', 2);

hold on

% Desenhar linha horizontal
plot([x_esq x_dir], [y_linha y_linha], 'k-', 'LineWidth', 2,'handlevisibility','off');

% Ajustes visuais
xlim([x_esq-1, x_dir+1]);
ylim([y_inf-1, y_sup+1]);
grid on;

hold off;