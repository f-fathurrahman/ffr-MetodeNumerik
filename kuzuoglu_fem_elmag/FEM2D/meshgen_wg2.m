function [conn,z,x,pmlbinr_node,pmlboutr_node,pmlr_node, ...
    pmlbinl_node,pmlboutl_node,pmll_node,scatb_node,btop_node, ...
    bbottom_node,scatin_node,scatin_elm,scatb_elm,dz1_node,dz2_node] ...
    = meshgen_wg2(obj)  
% Mesh generation function for "Scattering from an object within a parallel
% plate waveguide" using triangular elements
% Used in FEM2D_wg2.m
%
% INPUT:
% obj  : object structure (see the content of the structure in the main file)
% OUTPUT:
% conn : connectivity matrix (size: Mx3)
% z,x  : x and z coordinates of all nodes in the mesh (each size: Nx1)
% pmlbinr_node : Array containing the nodes on the inner PML boundary (right)
% pmlboutr_node : Array containing the nodes on the outer PML boundary (right)
% pmlr_node : Array containing the nodes inside the PML region (right)
% pmlbinl_node : Array containing the nodes on the inner PML boundary (left)
% pmlboutl_node : Array containing the nodes on the outer PML boundary (left)
% pmll_node : Array containing the nodes inside the PML region (left)
% scatb_node : Array containing the nodes on the surface of the scatterer
% btop_node : Array containing the nodes on the top waveguide wall
% bbottom_node : Array containing the nodes on the bottom waveguide wall
% scatin_node : Array containing the nodes inside the scatterer (empty for PEC scatterer)
% scatin_elm : Array containing the elements inside the scatterer (empty for PEC scatterer)
% scatb_elm : Array containing the elements on the surface of the scatterer
% dz1_node: Array containing the nodes on the measurement boundary on the left
% dz2_node: Array containing the nodes on the measurement boundary on the right
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
delh = obj.delh;             % element size
scat_type = obj.scat_type;   % 'pec' or 'diel' (PEC or dielectric object)
dx = obj.dx;                 % length of the waveguide along x
dz = obj.dz;                 % length of the waveguide along z
scat_shape = obj.scat_shape; % 1: circular, 2: rectangular
rc = obj.rc;                 % radius of circular object
lz = obj.lz;                 % half-length of rectangular object along z
lx = obj.lx;                 % half-length of rectangular object along x
pmldim = obj.pmldim;         % distance btw inner and outer PML boundaries
dz1 = obj.dz1;               % z1 place for reflection coefficient measurement
dz2 = obj.dz2;               % z2 place for transmission coefficient measurement

cx = obj.cx;                 % center of the object (along x)
cz = dz/2;                   % center of the object (along z)

% fit to delh incrementation
rc = round(rc/delh)*delh; 
lz = round(lz/delh)*delh;
lx = round(lx/delh)*delh;
cx = round(cx/delh)*delh; 
dx = round(dx/delh)*delh; 
dz = round(dz/delh)*delh;
pmldim = round(pmldim/delh)*delh; 
cz = round(cz/delh)*delh; 
cx = round(cx/delh)*delh; 
dz1 = round(dz1/delh)*delh;
dz2 = round(dz2/delh)*delh;

pmlinzr = dz; % right
pmloutzr = dz + pmldim; % right
pmlinzl = 0; % left
pmloutzl = 0 - pmldim; % left

% *************************************************************************
% Create coordinates 
% *************************************************************************
Nz = round((pmloutzr-pmloutzl)/delh)+1; 
Nx = round((dx-0)/delh)+1; 

zt = linspace(pmloutzl, pmloutzr, Nz);
xt = linspace(0, dx, Nx);

[za, xa] = meshgrid(zt, xt);
zat = za'; xat = xa';
z = zat(:); x = xat(:);
clear zat; clear xat;

if scat_shape == 1
    N = floor(pi/(2*asin(delh*0.5/rc)));
    alfa = cat(2, (pi/N)*(0:N), -fliplr((pi/N)*(1:N-1)));
    scatb_z = rc*cos(alfa)' + cz; % z-coordinates of the circle
    scatb_x = rc*sin(alfa)' + cx; % x-ccordinates of the circle
    
    threshold = delh*0.8;
    [discard,d] = knnsearch([z x],[scatb_z scatb_x],'k',10);
    ind = find(d < threshold);
    discard = unique(discard(ind));
    nodes0 = ones(1,length(z));
    nodes0(discard) = 0;
    nondiscard = find(nodes0 == 1);
    z = z(nondiscard);
    x = x(nondiscard);
    z = cat(1, z, scatb_z);
    x = cat(1, x, scatb_x);
    
    coords = [z x];
    coords = unique(coords,'rows');
    z = coords(:,1);
    x = coords(:,2);
    
    if strcmpi(scat_type, 'pec')    % PEC
        [IN, ON] = inpolygon(z, x, scatb_z, scatb_x);
        discard_node = find((IN == 1) & (ON == 0));
        
        nodes0 = ones(1,length(z));
        nodes0(discard_node) = 0;
        nondiscard_node = find(nodes0 ~= 0);
        z = z(nondiscard_node);
        x = x(nondiscard_node);
        z = [z;mean(scatb_z)];
        x = [x;mean(scatb_x)];
    end
    
else % scat_shape == 2
    scatb_node = find(((z < lz+cz+1e-8) & (z > lz+cz-1e-8) & ...
                       (x < lx+cx+1e-8) & (x > -lx+cx-1e-8)) | ...
                      ((z < -lz+cz+1e-8) & (z > -lz+cz-1e-8) & ...
                       (x < lx+cx+1e-8) & (x > -lx+cx-1e-8)) | ...
                      ((x < lx+cx+1e-8) & (x > lx+cx-1e-8) & ...
                       (z < lz+cz+1e-8) & (z > -lz+cz-1e-8)) | ...
                      ((x < -lx+cx+1e-8) & (x > -lx+cx-1e-8) & ...
                       (z < lz+cz+1e-8) & (z > -lz+cz-1e-8)));
    scatb_z = z(scatb_node);
    scatb_x = x(scatb_node);    
end

% *************************************************************************
% Create mesh and connectivity matrix
% *************************************************************************
conn = delaunay(z, x);    

N = length(z);
M = size(conn,1);

% *************************************************************************
% Find discarded elements and modify mesh if scattterer is PEC
% (make inner part of scatterer empty)
% *************************************************************************
if strcmpi(scat_type, 'pec')    % PEC
    zmid = mean(z(conn),2);
    xmid = mean(x(conn),2);
    
    if scat_shape == 1
       INPO = inpolygon(zmid, xmid, scatb_z, scatb_x);
    else
       INPO = inpolygon(zmid, xmid, [-lz lz lz -lz]+cz, [-lx -lx lx lx]+cx);        
    end
    discard_elm = find(INPO == 1);
    
    elms0 = ones(1,M);
    elms0(discard_elm) = 0;
    nondiscard_elm = find(elms0 ~= 0);
    
    conn = conn(nondiscard_elm,:);
    
    nodes = ones(1,N);
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

N = length(z);
M = size(conn,1);
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

zmid = mean(z(conn),2);
xmid = mean(x(conn),2);

% *************************************************************************
% Find  elements and nodes inside the scatterer if dielectric
% *************************************************************************
scatin_elm = [];  scatin_node = [];
if strcmpi(scat_type, 'diel') % diel
    if scat_shape == 1
       INPO = inpolygon(zmid, xmid, scatb_z, scatb_x);
    else
       INPO = inpolygon(zmid, xmid, [-lz lz lz -lz]+cz, [-lx -lx lx lx]+cx);        
    end    
    scatin_elm = find(INPO == 1);    
    
    if scat_shape == 1
       INPO = inpolygon(z, x, scatb_z, scatb_x);
    else
       INPO = inpolygon(z, x, [-lz-1e-8 lz+1e-8 lz+1e-8 -lz-1e-8]+cz, [-lx-1e-8 -lx-1e-8 lx+1e-8 lx+1e-8]+cx);        
    end     
    scatin_node = find(INPO == 1);
end

% *************************************************************************
% Find boundary nodes
% *************************************************************************
% nodes on the scatterer boundary
[scatb_node,~] = knnsearch([z x],[scatb_z scatb_x],'k',1);
scatb_node = sort(unique(scatb_node));

% right PML
pmlbinr_node = find((z > pmlinzr-1e-8) & (z < pmlinzr+1e-8));
pmlboutr_node = find((z > pmloutzr-1e-8) & (z < pmloutzr+1e-8));
pmlr_node = find((z > pmlinzr+1e-8) & (z < pmloutzr+1e-8));

% left PML
pmlbinl_node = find((z > pmlinzl-1e-8) & (z < pmlinzl+1e-8));
pmlboutl_node = find((z > pmloutzl-1e-8) & (z < pmloutzl+1e-8));
pmll_node = find((z < pmlinzl-1e-8) & (z > pmloutzl-1e-8));

% waveguide boundaries
btop_node = find((x > dx-1e-8) & (x < dx+1e-8));
bbottom_node = find((x > 0-1e-8) & (x < 0+1e-8));

% measurement boundaries
dz1_node = find((z > dz1-1e-8) & (z < dz1+1e-8));
dz2_node = find((z > dz2-1e-8) & (z < dz2+1e-8));

% *************************************************************************
% Find boundary elements on the scatterer if PEC
% *************************************************************************
scatb_elm = []; 
if strcmpi(scat_type, 'pec') % PEC
    nodesall = zeros(1, N);
    nodesall(scatb_node) = 1;
    
    ii = 1;
    for e = 1:M
        
        count = 0;
        for i = 1:3
            if (nodesall(conn(e,i)) == 1)
                count = count+1;
            end
        end
        
        if (count == 2)
            scatb_elm(ii) = e;
            ii = ii+1;
        end
    end
end

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
    
    plot(z(scatb_node), x(scatb_node),'r.', 'LineWidth',1.5)
    plot(z(pmlboutr_node), x(pmlboutr_node),'g.', 'LineWidth',1.5)
    plot(z(pmlbinr_node), x(pmlbinr_node),'m.', 'LineWidth',1.5)
    plot(z(pmlr_node), x(pmlr_node), 'y.', 'LineWidth',1.5)
    plot(z(pmlboutl_node), x(pmlboutl_node),'r.', 'LineWidth',1.5)
    plot(z(pmlbinl_node), x(pmlbinl_node),'b.', 'LineWidth',1.5)
    plot(z(pmll_node), x(pmll_node), 'c.', 'LineWidth',1.5)    
    plot(z(btop_node), x(btop_node), 'k.', 'LineWidth',1.5)  
    plot(z(bbottom_node), x(bbottom_node), 'k.', 'LineWidth',1.5)     
    plot(z(scatin_node), x(scatin_node), 'g.', 'LineWidth',1.5)
    plot(zmid(scatin_elm), xmid(scatin_elm), 'y.') 
    plot(zmid(scatb_elm), xmid(scatb_elm), 'r.') 
   
    plot(z(dz1_node), x(dz1_node),'m.', 'LineWidth',1.5)
    plot(z(dz2_node), x(dz2_node),'m.', 'LineWidth',1.5)
    
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
