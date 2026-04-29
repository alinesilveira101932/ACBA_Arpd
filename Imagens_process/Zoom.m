% Abrir figura
 openfig('Incertezas_pert.fig');
ax1 = gca;
% figure (2)
% Criar inset
ax2 = axes('Position', [0.43 0.53 0.4 0.25]);
copyobj(ax1.Children, ax2);

% Definir zoom (exemplo: em torno de x=2 a 4, y=0 a 2)
xlim(ax2, [0 10]);
ylim(ax2, [-0.1 0.5]);
grid(ax2, 'on'); box on; hold on;  
% plot([0 10],[0 0],'--r','linewidth',1,'handlevisibility','off')
% title(ax2, 'Zoom');

% Destacar área no gráfico principal
hold(ax1, 'on');
% rectangle(ax1, 'Position', [2, -1.5, 2.5, 2.5], 'EdgeColor', 'k', 'LineStyle',':','LineWidth', 1);
