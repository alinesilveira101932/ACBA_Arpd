clear, clc, close all

anom = 35.48;
% MPC azul
amin1 = 35.13;
amax1 = 53.22;
% MPC-RM magenta
amin2 = 5.32;
amax2 = 56.76;

bnom = 0.87;
bmin1 = -2.9;
bmax1 = 1.17;
bmin2 = -3.25;
bmax2 = 200;

cnom = 1.09*10^(-6);
cmin1 = 6.54*10^(-7);
cmax1 = 1.10*10^(-6);
cmin2= 6.54*10^(-7);
cmax2 = 2.62*10^(-4);

figure

subplot(1,3,1)
rectangleChatGPT(amin1, amax1, 1, 0.5, anom,[0.01,0.15,0.59])
hold on
rectangleChatGPT(amin2, amax2, 2, 0.5, anom, [0.72,0.27,1.00])
box on
xlim([0.5 2.5])
ylabel('a (s^{-2})','FontSize',11)
h = gca;
set(h,'XGrid','off','XTick',[])

subplot(1,3,2)
rectangleChatGPT(bmin1, bmax1, 1, 0.5, bnom,[0.01,0.15,0.59])
hold on
rectangleChatGPT(bmin2, bmax2, 2, 0.5, bnom, [0.72,0.27,1.00])
box on
xlim([0.5 2.5]); 
ylabel('b (s^{-1})','FontSize',11)
h = gca;
set(h,'XGrid','off','XTick',[])

subplot(1,3,3)
rectangleChatGPT(cmin1, cmax1, 1, 0.5, cnom,[0.01,0.15,0.59])
hold on
rectangleChatGPT(cmin2, cmax2, 2, 0.5, cnom, [0.72,0.27,1.00])
box on
xlim([0.5 2.5]); yscale('log')
ylabel('c (rad s^{-2} rpm^{-2})','FontSize',11)

hold on
% XTick: [0.5000 1 1.5000 2 2.5000]
plot([1 1],[90 90],'-','Color',[0.01,0.15,0.59],'LineWidth', 2)
plot([1 1],[90 90],'-','Color',[0.72,0.27,1.00],'LineWidth', 2)
plot([1 1],[90 90],'-k','LineWidth', 2)
 ylim([0 2.65*10^(-4)])
legend('MPC', 'MPC-RM', 'valor nominal','FontSize',12,'Location','southout','Orientation','horizontal')
h = gca;
set(h,'XGrid','off','XTick',[])