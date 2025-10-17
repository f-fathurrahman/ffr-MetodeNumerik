function [conn,z,x,pmlbin_node,pmlbout_node,pml_node,sourceb_node, ...
    btop_node,bbottom_node] = meshgen_wg1a(obj)
% Mesh generation function for "Example 1: Parallel plate waveguide with
% step discontinuity" using triangular elements
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
lz = obj.lz;      % length of the discontinuity along z
lx = obj.lx;      % length of the discontinuity along x
discontinuity_exist = obj.discontinuity_exist; % flag showing whether the discontinuity exists or not
pmldim = obj.pmldim;  % distance btw inner and outer PML boundaries

% fit to delh incrementation
lz = round(lz/delh)*delh; 
lx = round(lx/delh)*delh; 
dx = round(dx/delh)*delh; 
dz = round(dz/delh)*delh; 
pmldim = round(pmldim/delh)*delh;

pmlinz = dz;
pmloutz = dz + pmldim;

% *************************************************************************
% Create coordinates 
% *************************************************************************
Nz = round((pmloutz-0)/delh)+1; 
Nx = round((dx-0)/delh)+1; 

zt = linspace(0, pmloutz, Nz);
xt = linspace(0, dx, Nx);

[za, xa] = meshgrid(zt, xt);
zat = za'; xat = xa';
z = zat(:); x = xat(:);
clear zat; clear xat;

% *************************************************************************
% Create mesh and connectivity matrix
% *************************************************************************
conn = delaunay(z, x);    

% *************************************************************************
% Find discarded elements and modify mesh if discontinuity exists
% *************************************************************************
zmid = mean(z(conn),2);
xmid = mean(x(conn),2);
zc = mean(z);
zc = round(zc/delh)*delh; 

if discontinuity_exist
    dlx = round(((dx-lx)/2)/delh)*delh; 
    
    INPO = inpolygon(zmid, xmid, [zc-lz/2 zc+lz/2 zc+lz/2 zc-lz/2], [0 0 dlx dlx]);
    discard_elm1 = find(INPO == 1);
    INPO = inpolygon(zmid, xmid, [zc-lz/2 zc+lz/2 zc+lz/2 zc-lz/2], [lx+dlx lx+dlx dx dx]);
    discard_elm2 = find(INPO == 1);
        
    elms0 = ones(1,size(conn,1));
    elms0(discard_elm1) = 0;
    elms0(discard_elm2) = 0;
    nondiscard_elm = find(elms0 ~= 0);
    
    conn = conn(nondiscard_elm,:);
    
    nodes = ones(1,length(z));
    nodes(conn(:,1)) = 0;
    nodes(conn(:,2)) = 0;
    nodes(conn(:,3)) = 0;
    nondiscard_node = find(nodes == 0);
    discard_node = find(nodes ~= 0);
    
    z = z(nondiscard_node);
    x = x(nondiscard_node);
    
    % renumbering
    if ~isempty(discard_node)
        connd = conn(:);
        [connd, indx] = sort(connd);
        szc = size(connd);
        connd = diff(connd);
        connd(connd>0) = 1;
        connd = cumsum(connd)+1;
        connd = [1; connd];
        connd2 = zeros(szc);
        connd2(indx) = connd;
        clear connd;
        conn = reshape(connd2, size(conn));
        clear connd2;
    end
end

N = length(z);    % number of nodes
M = size(conn,1); % number of elements
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

% *************************************************************************
% Find boundary nodes
% *************************************************************************
sourceb_node = find((z > 0-1e-8) & (z < 0+1e-8)); % nodes on the boundary where the initial field is imposed

pmlbin_node = find((z > pmlinz-1e-8) & (z < pmlinz+1e-8));    % nodes on the inner PML boundary
pmlbout_node = find((z > pmloutz-1e-8) & (z < pmloutz+1e-8)); % nodes on the outer PML boundary
pml_node = find((z > pmlinz+1e-8) & (z < pmloutz+1e-8));      % nodes inside the PML region

btop_node = find((x > dx-1e-8) & (x < dx+1e-8));   % nodes on the waveguide top boundary
bbottom_node = find((x > 0-1e-8) & (x < 0+1e-8));  % nodes on the waveguide bottom boundary

if discontinuity_exist
    bbottom_node2 = find((z > zc-lz/2-1e-8) & (z < zc-lz/2+1e-8) & (x > 0-1e-8) & (x < dlx+1e-8));
    bbottom_node3 = find((z > zc+lz/2-1e-8) & (z < zc+lz/2+1e-8) & (x > 0-1e-8) & (x < dlx+1e-8));
    bbottom_node4 = find((x > dlx-1e-8) & (x < dlx+1e-8) & (z > zc-lz/2-1e-8) & (z < zc+lz/2+1e-8));
    
    btop_node2 = find((z > zc-lz/2-1e-8) & (z < zc-lz/2+1e-8) & (x > lx+dlx-1e-8) & (x < dx+1e-8));
    btop_node3 = find((z > zc+lz/2-1e-8) & (z < zc+lz/2+1e-8) & (x > lx+dlx-1e-8) & (x < dx+1e-8));
    btop_node4 = find((x > lx+dlx-1e-8) & (x < lx+dlx+1e-8) & (z > zc-lz/2-1e-8) & (z < zc+lz/2+1e-8));
    
    btop_node = [btop_node; btop_node2; btop_node3; btop_node4];
    bbottom_node = [bbottom_node; bbottom_node2; bbottom_node3; bbottom_node4];
end

nodes0 = zeros(N,1); 
nodes0(btop_node) = 1;
nodes0(sourceb_node) = 0;
btop_node = find(nodes0 == 1);

nodes0 = zeros(N,1); 
nodes0(bbottom_node) = 1;
nodes0(sourceb_node) = 0;
bbottom_node = find(nodes0 == 1);

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
    
    plot(z(sourceb_node), x(sourceb_node),'r.', 'LineWidth',1.5)
    plot(z(pmlbout_node), x(pmlbout_node),'g.', 'LineWidth',1.5)
    plot(z(pmlbin_node), x(pmlbin_node),'m.', 'LineWidth',1.5)
    plot(z(pml_node), x(pml_node), 'y.', 'LineWidth',1.5)
    plot(z(btop_node), x(btop_node), 'k.', 'LineWidth',1.5)  
    plot(z(bbottom_node), x(bbottom_node), 'k.', 'LineWidth',1.5)      
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
