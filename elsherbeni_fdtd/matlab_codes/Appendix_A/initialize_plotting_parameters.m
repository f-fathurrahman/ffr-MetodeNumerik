% subroutine used to initialize 1D plot

Ez_positions = [0:nx]*dx;
Hy_positions = ([0:nx-1]+0.5)*dx;
v = [0 -0.1 -0.1; 0 -0.1 0.1; 0 0.1 0.1; 0 0.1 -0.1; ...
     1 -0.1 -0.1; 1 -0.1 0.1; 1 0.1 0.1; 1 0.1 -0.1];
f = [1 2 3 4; 5 6 7 8];
axis([0 1 -0.2 0.2 -0.2 0.2]);
lez = line(Ez_positions,Ez*0,Ez,'Color','b','LineWidth',1.5);
lhy = line(Hy_positions,377*Hy,Hy*0,'Color','r', ...
    'LineWidth',1.5,'linestyle','-.');
set(gca,'fontsize',12,'FontWeight','bold');
axis square;
legend('E_{z}', 'H_{y} \times 377','Location','NorthEast');
xlabel('x [m]');
ylabel('[A/m]');
zlabel('[V/m]');
grid on;
p = patch('vertices', v, 'faces', f, 'facecolor', 'g', 'facealpha',0.2);
text(0,1,1.1,'PEC','horizontalalignment','center','fontweight','bold');
text(1,1,1.1,'PEC','horizontalalignment','center','fontweight','bold');
