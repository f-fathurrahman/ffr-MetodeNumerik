%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%    An Introduction to Scientific Computing          %%%%%%%
%%%%%%%    I. Danaila, P. Joly, S. M. Kaber & M. Postel     %%%%%%%
%%%%%%%                 Springer, 2023                      %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Visualization of iso-contours                 %
%   of the passive scalar (or tracer)             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NSE_F_visu_sca(x,y,sca,niso,Tstep,time)

global dx dy Lx Ly 
global im ip jp jm ic jc

fs =20; % font size
lw =2;  % line width
mk =10; % markersize

[xx,yy]=meshgrid(x,y);xx=xx';yy=yy';

colormap jet
pcolor(xx,yy,sca);
axis([0,Lx,0,Ly]);axis equal;
shading('interp');drawnow
xlabel('x');ylabel('y');
title(['Scalar   t=' num2str(time) ', Tstep=' num2str(Tstep)],'FontSize',fs);
set(gca,'FontSize',fs,'Xtick',0:Lx/4:Lx,'Ytick',0:Ly/4:Ly);
