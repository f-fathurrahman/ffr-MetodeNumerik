function [conn,z,x,pmlbin_node,pmlbout_node,pml_node,sourceb_node, ...
    btop_node,bbottom_node] = meshgen_wg1b(obj)
% Mesh generation function for "Example 2: Parallel plate waveguide mounted
% on ground" using triangular elements
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
lx = obj.lx;      % length of the ground along x
pmldim = obj.pmldim;  % distance btw inner and outer PML boundaries
fsdim = obj.fsdim;    % length of the free space region along z (btw ground and inner PML boundary)

% fit to delh incrementation
lx = round(lx/delh)*delh; 
dx = round(dx/delh)*delh; 
dz = round(dz/delh)*delh; 
pmldim = round(pmldim/delh)*delh;
fsdim = round(fsdim/delh)*delh;

% create critical points
xp(1) = dx; zp(1) = 0;
xp(2) = xp(1); zp(2) = dz;
xp(3) = xp(2)+lx+pmldim; zp(3) = zp(2);
xp(4) = xp(3); zp(4) = zp(3)+fsdim+pmldim;
xp(5) = -(lx+pmldim); zp(5) = zp(4);
xp(6) = xp(5); zp(6) = zp(2);
xp(7) = 0; zp(7) = zp(6);
xp(8) = xp(7); zp(8) = 0;

% *************************************************************************
% Create coordinates 
% *************************************************************************
Nz = round(zp(5)/delh)+1; 
Nx = round((xp(3)-xp(6))/delh)+1; 

zt = linspace(0, zp(5), Nz);
xt = linspace(xp(6), xp(3), Nx);

[za, xa] = meshgrid(zt, xt);
zat = za'; xat = xa';
z = zat(:); x = xat(:);
clear zat; clear xat;

% create boundary
zpv = [zp zp(1)];
xpv = [xp xp(1)];

svdist = sqrt((xpv(2:end)-xpv(1:end-1)).^2+(zpv(2:end)-zpv(1:end-1)).^2);
nd = round(svdist/delh)+1;
NV = length(svdist);

bz = []; bx = [];
for i = 1:NV
    tempz = linspace(zpv(i), zpv(i+1), nd(i));
    tempx = linspace(xpv(i), xpv(i+1), nd(i));
    bz = [bz; tempz'];
    bx = [bx; tempx'];
end

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

% *************************************************************************
% Find discarded elements and modify mesh if discontinuity exists
% *************************************************************************
zmid = mean(z(conn),2);
xmid = mean(x(conn),2);

INPO = inpolygon(zmid, xmid, zpv,xpv);
discard_elm1 = find(INPO == 0);

elms0 = ones(1,size(conn,1));
elms0(discard_elm1) = 0;
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

N = length(z);    % number of nodes
M = size(conn,1); % number of elements
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

% *************************************************************************
% Find boundary nodes
% *************************************************************************
% nodes on the boundary where the initial field is imposed
sourceb_node = find((z > 0-1e-8) & (z < 0+1e-8)); 

% nodes on the inner PML boundary
pmlbin_node = find(((x > xpv(3)-pmldim-1e-8) & (x < xpv(4)-pmldim+1e-8) & (z > zpv(3)-1e-8) & (z < zpv(4)-pmldim+1e-8)) | ...
    ((x > xpv(5)+pmldim-1e-8) & (x < xpv(4)-pmldim+1e-8) & (z > zpv(5)-pmldim-1e-8) & (z < zpv(4)-pmldim+1e-8)) | ...
    ((x > xpv(5)+pmldim-1e-8) & (x < xpv(6)+pmldim+1e-8) & (z > zpv(6)-1e-8) & (z < zpv(5)-pmldim+1e-8)));    

% nodes on the outer PML boundary
pmlbout_node = find(((x > xpv(3)-1e-8) & (x < xpv(4)+1e-8) & (z > zpv(3)-1e-8) & (z < zpv(4)+1e-8)) | ...
    ((x > xpv(5)-1e-8) & (x < xpv(4)+1e-8) & (z > zpv(5)-1e-8) & (z < zpv(4)+1e-8)) | ...
    ((x > xpv(5)-1e-8) & (x < xpv(6)+1e-8) & (z > zpv(6)-1e-8) & (z < zpv(5)+1e-8)));    

% nodes inside the PML region
pml_node = find(((x > xpv(3)-pmldim+1e-8) & (x < xpv(3)+1e-8)) | ...
    ((x > xpv(6)+1e-8) & (x < xpv(6)+pmldim-1e-8)) | ...
    ((z > zpv(5)-pmldim+1e-8) & (z < zpv(5)+1e-8)));   

% nodes on the waveguide top boundary
btop_node = find(((x > xpv(1)-1e-8) & (x < xpv(2)+1e-8) & (z > zpv(1)-1e-8) & (z < zpv(2)+1e-8)) | ...
                 ((x > xpv(2)-1e-8) & (x < xpv(3)+1e-8) & (z > zpv(2)-1e-8) & (z < zpv(3)+1e-8)));   

% nodes on the waveguide bottom boundary
bbottom_node = find(((x > xpv(8)-1e-8) & (x < xpv(7)+1e-8) & (z > zpv(8)-1e-8) & (z < zpv(7)+1e-8)) | ...
                 ((x > xpv(6)-1e-8) & (x < xpv(7)+1e-8) & (z > zpv(6)-1e-8) & (z < zpv(7)+1e-8)));   

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
