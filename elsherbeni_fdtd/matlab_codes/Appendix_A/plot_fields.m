% subroutine used to plot 1D transient fields

delete(lez);
delete(lhy);
lez = line(Ez_positions,Ez*0,Ez,'Color','b','LineWidth',1.5);
lhy = line(Hy_positions,377*Hy,Hy*0,'Color','r', ...
    'LineWidth',1.5,'linestyle','-.');
ts = num2str(time_step);
ti = num2str(dt*time_step*1e9);
title(['time step = ' ts ', time = ' ti ' ns']);
drawnow;
