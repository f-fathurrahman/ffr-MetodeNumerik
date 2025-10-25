function tsolmovie(S,tt,xy,x,y,fig,cmap,axis_label)
%TSOLMOVIE plots movie for Q1 solution data
%   tsolmovie(U,time,xy,ev,fig,cmap,axis_label);
%   input
%          U           nodal solution vector
%          time        time vector
%          xy          nodal coordinate vector
%          x, y        plotting grid vectors
%          fig         figure number
%          cmap        colormap (optional; default is 'winter')
%          axis_label  optional axis values for plotting
%
%   IFISS function: DJS; 2 June 2018.
% Copyright (c) 2011 D.J. Silvester
if nargin<7, cmap='winter'; end
if nargin<8,
axis([min(xy(:,1)),max(xy(:,1)),min(xy(:,2)),max(xy(:,2)), ...
            0,max(S(:,1))]); end
fprintf('plotting solution movie... ')
av=axis;
[ndof,nsnaps]=size(S);
[X,Y]=meshgrid(x,y);
figure(fig)
colormap(cmap)
for snap=2:nsnaps
sol=S(:,snap); ttk=tt(snap);
% interpolate to a cartesian product mesh
xysol = griddata(xy(:,1),xy(:,2),sol,X,Y);
subplot(122),contour(X,Y,xysol,15),  axis('square'), 
if all([min(x),max(x),min(y),max(y)] == [-1,1,-1,1]), squarex, end
title(['time is ',num2str(ttk,'%7.4f')],'FontSize',12),
subplot(121),mesh(X,Y,xysol),axis('square'), axis(av),
view(330,30)
title(['time is ',num2str(ttk,'%7.4f')],'FontSize',12)
pause(0.1)
end
fprintf('done\n')
return
