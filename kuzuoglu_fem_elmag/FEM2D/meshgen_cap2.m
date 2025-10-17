function [conn,x,y,xmid,ymid,btop_node,bbottom_node,bleft_node,bright_node,d_elm] = meshgen_cap2(Lx,Ly,dx,dy,delx,dely)
% Mesh generation function for "Parallel-Plate Capacitor Problem with Dielectric Inhomogeneity" 
% with triangular elements
% Used in FEM2D_cap2.m
%
% INPUT:
% Lx : length of the computational domain along x axis (meter)
% Ly : length of the computational domain along y axis (meter)
% dx : length of the dielectric block along x axis (meter)
% dy : length of the dielectric block along y axis (meter)
% delx : increment along x axis (meter)
% dely : increment along y axis (meter)
% OUTPUT:
% conn : connectivity matrix (size: Mx3)
% x,y  : x and y coordinates of all nodes in the mesh (each size: Nx1)
% xmid,ymid : x and y midpoints of all elements (each size: Mx1)
% btop_node : Array containing the nodes on the top plate
% bbottom_node : Array containing the nodes on the bottom plate
% bleft_node : Array containing the nodes on the left outer boundary
% bright_node : Array containing the nodes on the right outer boundary
% d_elm : Array containing the elements within the dielectric block
%
% Note: M is the number of elements, N is the number of nodes
%
%
% Copyright (C) 2018, Ozlem Ozgun, Mustafa Kuzuoglu
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%    
% Contact:
% Dr. Ozlem Ozgun
% email: ozgunozlem@gmail.com
% web: http://www.ee.hacettepe.edu.tr/~ozlem/

% fit to delh incrementation
Lx = round(Lx/delx)*delx; 
Ly = round(Ly/dely)*dely; 
dx = round(dx/delx)*delx; 
dy = round(dy/dely)*dely; 

% *************************************************************************
% Create coordinates 
% *************************************************************************
Nx = round(Lx/delx)+1; 
Ny = round(Ly/dely)+1; 

xt = linspace(-Lx/2, Lx/2, Nx);
yt = linspace(-Ly/2, Ly/2, Ny);   

[xa, ya] = meshgrid(xt, yt);
xa = xa'; ya = ya';
x = xa(:); y = ya(:);

% *************************************************************************
% Create mesh and connectivity matrix
% *************************************************************************
conn = delaunay(x, y);

N = length(x);
M = size(conn,1);
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

% *************************************************************************
% Find nodes and elements
% *************************************************************************
% Nodes on the outer boundary
btop_node = find((y < Ly/2+1e-8) & (y > Ly/2-1e-8) & (x < Lx/2+1e-8) & (x > -Lx/2-1e-8));
bbottom_node = find((y < -Ly/2+1e-8) & (y > -Ly/2-1e-8) & (x < Lx/2+1e-8) & (x > -Lx/2-1e-8));
bleft_node = find((x < -Lx/2+1e-8) & (x > -Lx/2-1e-8) & (y < Ly/2+1e-8) & (y > -Ly/2-1e-8));
bright_node = find((x < Lx/2+1e-8) & (x > Lx/2-1e-8) & (y < Ly/2+1e-8) & (y > -Ly/2-1e-8));

% Coordinates of the midpoints of the elements
xmid = mean(x(conn),2);
ymid = mean(y(conn),2);

% Elements inside the dielectric block
d_elm = find((ymid < dy/2) & (ymid > -dy/2) & (xmid > -dx/2) & (xmid < dx/2));

% *************************************************************************
% Plot results
% *************************************************************************
plot_flag = 1; % flag showing whether the mesh will be plotted or not
if plot_flag
    figure;
    triplot(conn, x, y)
    axis equal tight;
    set(gcf,'Color',[1 1 1]);
    xlabel('x (m)');  ylabel('y (m)');
    hold on
    plot(x(btop_node), y(btop_node),'k.','linewidth',2)
    plot(x(bbottom_node), y(bbottom_node),'r.','linewidth',2)
    plot(x(bleft_node), y(bleft_node),'g.','linewidth',2)
    plot(x(bright_node), y(bright_node),'m.','linewidth',2)    
    plot(xmid(d_elm), ymid(d_elm),'y.')
end

% *************************************************************************
% Calculate element quality
% *************************************************************************
qual = triangle_quality2(x(conn),y(conn));

figure; set(gcf,'Color',[1 1 1]);
plot(1:M, qual, 'k'); 
title('Element Quality');
xlabel('Element No');
ylabel('Quality');
hold on
ys = 0.25*ones(1,M);
plot(ys, 'r','LineWidth',1.5)
axis([1 M 0 1]);
