function [conn,z,x,pmlbin_node,pmlbout_node,pml_node,sourceb_node, ...
    btop_node,bbottom_node] = meshgen_wg1d(obj)
% Mesh generation function for "Example 4: 2-D Horn" using triangular elements
% Used in FEM2D_wg1.m
%
% INPUT:
% obj  : object structure (see the content of the structure in the main file)
% OUTPUT:
% conn : connectivity matrix (size: Mx3)
% z,x  : x and z coordinates of all nodes in the mesh (each size: Nx1)
% pmlbin_node : Array containing the nodes on the inner PML boundary
% pmlbout_node : Array containing the nodes on the outer PML boundary
% pml_node : Array containing the nodes inside the PML region
% sourceb_node: Array containing the nodes on the boundary where the initial field is imposed
% btop_node : Array containing the nodes on the top waveguide wall
% bbottom_node : Array containing the nodes on the bottom waveguide wall
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

% Get input parameters
delh = obj.delh;  % element size
dz = obj.dz;      % length of the waveguide along z
dx = obj.dx;      % length of the waveguide along x
az = obj.az;      % length of the tilted edge of the horn along z
ax = obj.ax;      % length of the tilted edge of the horn along x
pmldim = obj.pmldim;  % distance btw inner and outer PML boundaries
fsdim = obj.fsdim;    % length of the free space region along z (btw ground and inner PML boundary)

% fit to delh incrementation
dx = round(dx/delh)*delh; 
dz = round(dz/delh)*delh; 
ax = round(ax/delh)*delh; 
az = round(az/delh)*delh; 
pmldim = round(pmldim/delh)*delh;
fsdim = round(fsdim/delh)*delh;

% create critical points
xp(1) = dx; zp(1) = 0;
xp(2) = xp(1); zp(2) = dz;
xp(3) = xp(2)+ax; zp(3) = zp(2)+az;

xp(4) = 0; zp(4) = 0;
xp(5) = 0; zp(5) = zp(2);
xp(6) = -ax; zp(6) = zp(3);

% *************************************************************************
% Create coordinates 
% *************************************************************************
Nz = round((dz+az+fsdim+pmldim)/delh)+1; 
Nx = round((dx+2*ax+2*fsdim+2*pmldim)/delh)+1; 

zt = linspace(0, dz+az+fsdim+pmldim, Nz);
xt = linspace(-(ax+fsdim+pmldim), (dx+ax+fsdim+pmldim), Nx);

[za, xa] = meshgrid(zt, xt);
zat = za'; xat = xa';
z = zat(:); x = xat(:);
clear zat; clear xat;

% create boundary
svdist = sqrt((xp(2:end)-xp(1:end-1)).^2+(zp(2:end)-zp(1:end-1)).^2);
nd = round(svdist/delh)+1;
NV = length(svdist);

btop_z = []; btop_x = []; bbottom_z = []; bbottom_x = [];
for i = 1:NV
    tempz = linspace(zp(i), zp(i+1), nd(i));
    tempx = linspace(xp(i), xp(i+1), nd(i));
    
    if i <= 2
       btop_z = [btop_z; tempz'];
       btop_x = [btop_x; tempx'];
    elseif (i >=  4)
       bbottom_z = [bbottom_z; tempz'];
       bbottom_x = [bbottom_x; tempx'];        
    end        
end
bz = [btop_z; bbottom_z];
bx = [btop_x; bbottom_x];

threshold = delh*0.8;
[discard,d] = knnsearch([z x],[bz bx],'k',10);
ind = find(d < threshold);
discard = unique(discard(ind));
nodes0 = ones(1,length(x));
nodes0(discard) = 0;
nondiscard = find(nodes0 == 1);
z = z(nondiscard);
x = x(nondiscard);
z = cat(1, z, bz);
x = cat(1, x, bx);
  
coords = [z x];
coords = unique(coords,'rows');
z = coords(:,1);
x = coords(:,2);
    
% *************************************************************************
% Create mesh and connectivity matrix
% *************************************************************************
conn = delaunay(z, x);    

N = length(z);    % number of nodes
M = size(conn,1); % number of elements
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

% *************************************************************************
% Find boundary nodes
% *************************************************************************
% nodes on the boundary where the initial field is imposed
sourceb_node = find((z > pmldim+delh-1e-8) & (z < pmldim+delh+1e-8) & (x > 0-1e-8) & (x < dx+1e-8)); 

% nodes on the waveguide top boundary
[btop_node,~] = knnsearch([z x],[btop_z btop_x],'k',1);
btop_node = sort(unique(btop_node));

% nodes on the waveguide bottom boundary
[bbottom_node,~] = knnsearch([z x],[bbottom_z bbottom_x],'k',1);
bbottom_node = sort(unique(bbottom_node));

pmlout_minx = -(ax+fsdim+pmldim);
pmlout_maxx = (dx+ax+fsdim+pmldim);
pmlout_minz = 0;
pmlout_maxz = dz+az+fsdim+pmldim;

pmlin_minx = pmlout_minx + pmldim;
pmlin_maxx = pmlout_maxx - pmldim;
pmlin_minz = pmlout_minz + pmldim;
pmlin_maxz = pmlout_maxz - pmldim;

% nodes on the inner PML boundary
pmlbin_node = find(((x < pmlin_maxx+1e-8) & (x > pmlin_maxx-1e-8) & ...
                    (z < pmlin_maxz+1e-8) & (z > pmlin_minz-1e-8)) | ...
                   ((x < pmlin_minx+1e-8) & (x > pmlin_minx-1e-8) & ...
                    (z < pmlin_maxz+1e-8) & (z > pmlin_minz-1e-8)) | ...
                   ((z < pmlin_maxz+1e-8) & (z > pmlin_maxz-1e-8) & ...
                    (x < pmlin_maxx+1e-8) & (x > pmlin_minx-1e-8)) | ...
                   ((z < pmlin_minz+1e-8) & (z > pmlin_minz-1e-8) & ...
                    (x < pmlin_maxx+1e-8) & (x > pmlin_minx-1e-8)));

% nodes on the outer PML boundary               
pmlbout_node = find(((x < pmlout_maxx+1e-8) & (x > pmlout_maxx-1e-8) & ...
                    (z < pmlout_maxz+1e-8) & (z > pmlout_minz-1e-8)) | ...
                   ((x < pmlout_minx+1e-8) & (x > pmlout_minx-1e-8) & ...
                    (z < pmlout_maxz+1e-8) & (z > pmlout_minz-1e-8)) | ...
                   ((z < pmlout_maxz+1e-8) & (z > pmlout_maxz-1e-8) & ...
                    (x < pmlout_maxx+1e-8) & (x > pmlout_minx-1e-8)) | ...
                   ((z < pmlout_minz+1e-8) & (z > pmlout_minz-1e-8) & ...
                    (x < pmlout_maxx+1e-8) & (x > pmlout_minx-1e-8)));

% nodes inside the PML region
INPO = inpolygon(z, x, [pmlin_minz pmlin_maxz pmlin_maxz pmlin_minz], [pmlin_maxx pmlin_maxx pmlin_minx pmlin_minx]);
temp_node = find(INPO == 1);      
nodes0 = ones(1,N); 
nodes0(temp_node) = 0;
pml_node = find(nodes0 == 1);

% *************************************************************************
% Plot results
% *************************************************************************
plot_flag = 1; % flag showing whether the mesh will be plotted or not
if plot_flag
    figure; hold on
    triplot(conn, z, x)
    axis equal tight;
    set(gcf,'Color',[1 1 1]);
    xlabel('z (m)');  ylabel('x (m)');
    
    plot(z(pmlbout_node), x(pmlbout_node),'g.', 'LineWidth',1.5)
    plot(z(pmlbin_node), x(pmlbin_node),'m.', 'LineWidth',1.5)
    plot(z(pml_node), x(pml_node), 'y.', 'LineWidth',1.5)
    plot(z(btop_node), x(btop_node), 'k.', 'LineWidth',1.5)  
    plot(z(bbottom_node), x(bbottom_node), 'k.', 'LineWidth',1.5)  
    plot(z(sourceb_node), x(sourceb_node),'r.', 'LineWidth',1.5)    
end

% *************************************************************************
% Calculate element quality
% *************************************************************************
qual = triangle_quality2(z(conn),x(conn));

figure; set(gcf,'Color',[1 1 1]);
plot(1:M, qual, 'k'); 
title('Element Quality');
xlabel('Element No');
ylabel('Quality');
hold on
ys = 0.25*ones(1,M);
plot(ys, 'r','LineWidth',1.5)
axis([1 M 0 1]);
