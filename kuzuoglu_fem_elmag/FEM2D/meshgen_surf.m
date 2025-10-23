function [conn,x,y,pmlbin_node,pmlbout_node,pml_node,...
      scatb_node,scatb_elm,huygb_node,huygb_elm,wind] ...
      = meshgen_surf(obj,freq)
% Mesh generation function for "Scattering from a Rough Surface"
% using triangular elements
% Used in FEM2D_surf.m
%
% INPUT:
% obj  : object structure (see the content of the structure in the main file)
% freq : frequency (Hz)
% OUTPUT:
% conn : connectivity matrix (size: Mx3)
% x,y  : x and y coordinates of all nodes in the mesh (each size: Nx1)
% pmlbin_node : Array containing the nodes on the inner PML boundary
% pmlbout_node : Array containing the nodes on the outer PML boundary
% pml_node : Array containing the nodes inside the PML region
% scatb_node : Array containing the nodes on the surface
% scatb_elm : Array containing the elements on the surface 
% huygb_node : Array containing the nodes on the surface of the Huygens boundary
% huygb_elm : Array containing the elements on the surface of the Huygens boundary
% wind: Window function for tapering the incident field on the surface
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

c0 = 3*1e8;         % m/sec, velocity of light in free space
lambda0 = c0/freq;  % meter, wavelength

delh = obj.delh;   % element size
Len = obj.Len;     % length of the rough surface
xdiscard = obj.xdiscard; % discarded region from the surface
surface_model = obj.surface_model; % 'flat', 'sea' or 'ground'
v = obj.v;   % wind velocity (m/s) (used for sea surface)
h = obj.h;   % rms height (used for ground surface)
lc = obj.lc; % correlation length (used for ground surface)

pmldim = 2.0*lambda0; % thickness of the PML region
fsdim = 2.0*lambda0;  % thickness of the free space region

% fit to delh incrementation
fsdim = round(fsdim/delh)*delh; 
pmldim = round(pmldim/delh)*delh; 

% generate the surface
[xsurf, ysurf] = rough_surface(delh, pmldim, Len, surface_model, v, h, lc);

ytop = round(max(ysurf)/delh)*delh + fsdim;

qual_thresh = 0.25;

pmlin_maxx = Len;
pmlin_maxy = ytop;
pmlin_minx = 0;
pmlin_miny = floor(min(ysurf)/delh)*delh; 

pmlout_maxx = pmlin_maxx + pmldim;
pmlout_maxy = pmlin_maxy + pmldim;
pmlout_minx = pmlin_minx - pmldim;
pmlout_miny = pmlin_miny;
          
disp('Mesh and Nodal Connectivity Matrix are being created ...');
 
Nx = round((pmlout_maxx-pmlout_minx)/delh)+1; 
Ny = round((pmlout_maxy-pmlout_miny)/delh)+1; 

xt = linspace(pmlout_minx, pmlout_maxx, Nx);
yt = linspace(pmlout_miny, pmlout_maxy, Ny);

[xa, ya] = meshgrid(xt, yt);
xat = xa'; yat = ya';
x = xat(:)'; y = yat(:)';
clear xat; clear yat;


%%
threshold = delh * 0.8;

if ~isempty(find(ysurf~=0))
    [discard,d] = knnsearch([x' y'],[xsurf' ysurf'],'k',40);
    ind = find(d < threshold);
    discard = unique(discard(ind));
    
    [discard1,d] = knnsearch([x' y'],[xsurf' (ysurf+delh)'],'k',40);
    ind = find(d < threshold);
    discard1 = unique(discard1(ind));
    
    [discard3,d] = knnsearch([x' y'],[xsurf' (ysurf+2*delh)'],'k',40);
    ind = find(d < threshold);
    discard3 = unique(discard3(ind));
    
    
    tempx = cat(2, xsurf, [pmlout_maxx pmlout_maxx pmlout_minx pmlout_minx]);
    tempy = cat(2, ysurf, [ysurf(end) min(y) min(y) ysurf(1)]);
    INPO = inpolygon(x, y, tempx, tempy);
    discard2 = find(INPO == 1);
    
    nodes0 = 1:length(x);
    nodes0(discard) = 0;
    nodes0(discard2) = 0;
    nodes0(discard1) = 0;
    nodes0(discard3) = 0;
    nondiscard = find(nodes0 ~= 0);
    
    x = x(nondiscard);
    y = y(nondiscard);
    
    x = cat(2, x, xsurf, xsurf, xsurf);
    y = cat(2, y, ysurf, ysurf+delh, ysurf+2*delh);
    
    coords = [x.' y.'];
    coords = unique(coords,'rows');
    x = coords(:,1).';
    y = coords(:,2).';
    
    clear coords;
    clear discard; clear discard1; clear discard2; clear discard3;
end

% *************************************************************************
% Create nodal connectivity matrix
% *************************************************************************
conn = delaunay(x, y);

% Find the center-of-mass points of elements 
xmid = mean(x(conn),2);
ymid = mean(y(conn),2);

% Find discarded elements and modify mesh 
discard_elm = [];
if ~isempty(find(ysurf~=0))
    tempx = cat(2, xsurf, [pmlout_maxx pmlout_maxx pmlout_minx pmlout_minx]);
    tempy = cat(2, ysurf, [ysurf(end) min(y) min(y) ysurf(1)]);
    INPO = inpolygon(xmid, ymid, tempx, tempy);
    discard_elm = find(INPO == 1);
end

% calculate element quality
qual = triangle_quality2(x(conn),y(conn));
        
discard_elm2 = find(qual < qual_thresh);
      
elms0 = 1:size(conn,1); 
elms0(discard_elm) = 0;
elms0(discard_elm2) = 0;
nondiscard_elm = find(elms0 ~= 0);

conn = conn(nondiscard_elm,:);

nodes = 1:length(x);
nodes(conn(:,1)) = 0;
nodes(conn(:,2)) = 0;
nodes(conn(:,3)) = 0;
nondiscard_node = find(nodes == 0);
discard_node = find(nodes ~= 0);

x = x(nondiscard_node);
y = y(nondiscard_node);

if ~isempty(discard_node)
    conncombined = conn(:);
    [connsorted, indx] = sort(conncombined);
    connd = diff(connsorted);
    connd(connd>0) = 1;
    conndc = cumsum(connd)+1;    
    connsorted2 = [1; conndc];    
    conncombined2 = zeros(size(conncombined));
    conncombined2(indx) = connsorted2;
    conn = reshape(conncombined2, size(conn));
end

clear elms0; clear qual;
clear nondiscard_elm;
clear discard_elm2; clear discard_elm3; 
clear nondiscard_node; clear discard_node;
 
N = length(y);
M = length(conn);
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);
disp(' ');


% *************************************************************************
% Find the center-of-mass points of elements
% *************************************************************************
xmid = mean(x(conn),2);
ymid = mean(y(conn),2);

% *************************************************************************
% Find scatterer boundary nodes
% *************************************************************************
pmlbin_node = find(((x < pmlin_maxx+1e-8) & (x > pmlin_maxx-1e-8) & ...
                    (y < pmlin_maxy+1e-8) & (y > pmlin_miny-1e-8)) | ...
                   ((x < pmlin_minx+1e-8) & (x > pmlin_minx-1e-8) & ...
                    (y < pmlin_maxy+1e-8) & (y > pmlin_miny-1e-8)) | ...
                   ((y < pmlin_maxy+1e-8) & (y > pmlin_maxy-1e-8) & ...
                    (x < pmlin_maxx+1e-8) & (x > pmlin_minx-1e-8)));
              
% *************************************************************************
% Find VOI boundary nodes and coordinates
% *************************************************************************
pmlbout_node = find(((x < pmlout_maxx+1e-8) & (x > pmlout_maxx-1e-8) & ...
                     (y < pmlout_maxy+1e-8) & (y > pmlout_miny-1e-8)) | ...
                    ((x < pmlout_minx+1e-8) & (x > pmlout_minx-1e-8) & ...
                     (y < pmlout_maxy+1e-8) & (y > pmlout_miny-1e-8)) | ...
                    ((y < pmlout_maxy+1e-8) & (y > pmlout_maxy-1e-8) & ...
                     (x < pmlout_maxx+1e-8) & (x > pmlout_minx-1e-8)));
                 
% *************************************************************************
% Find PML inside nodes
% *************************************************************************
pmlinside_node = find(((x < (pmlin_maxx-1e-8)) & (x > (pmlin_minx+1e-8))) & ...
                      ((y < (pmlin_maxy-1e-8)) & (y > (pmlin_miny+1e-8))));

pmlinsideb_node = cat(2, pmlinside_node, pmlbin_node);

ind = find((x < Len-1e-6) & (x > 0+1e-6) & (y == pmlin_miny));

nodes0 = 1:length(x); 
nodes0(pmlinsideb_node) = 0;
nodes0(ind) = 0;
pml_node = find(nodes0 ~= 0);
clear nodes0;


% *************************************************************************
% Find boundary nodes and inside elms
% *************************************************************************
[surf_node,~] = knnsearch([x' y'],[xsurf' ysurf'],'k',1);
surf_node = sort(unique(surf_node));
    
scatb_node = surf_node';
    
[val, ind] = min(abs(y-(max(ysurf)+4*delh)));
huygb_node = find((y < y(ind)+1e-8) & (y > y(ind)-1e-8));

xdiscard2 = xdiscard;
ind = find(~((x(huygb_node) < Len-xdiscard2-1e-6) & (x(huygb_node) > 0+xdiscard2+1e-6)));
huygb_node(ind) = [];

   
% *************************************************************************
% Find surface elements and nodes (needed for impedance surface)
% ************************************************************************* 
nodesall = zeros(1, N);
nodesall(huygb_node) = 1;

nc = nodesall(conn);
ncsum = sum(nc,2);
huygb_elm = find(ncsum == 2);
 
elm = zeros(1,M);
elm(huygb_elm) = 1; 
ind = find(ymid < y(huygb_node(1)));
elm(ind) = 0;
huygb_elm = find(elm==1); 


nodesall = zeros(1, N);    
nodesall(scatb_node) = 1;

nc = nodesall(conn);
ncsum = sum(nc,2);
scatb_elm = find(ncsum == 2);
scatb_elm = scatb_elm';

% ind = find(((xmid(scatb_elm) < Len-xdiscard2-1e-6) & (xmid(scatb_elm) > 0+xdiscard2+1e-6)));
% scatb_elm2 = scatb_elm(ind); 


%% ************************************************************************
% Window function
% *************************************************************************
hx = xdiscard;
tempx = real(x(surf_node)); tempy = real(y(surf_node));
xleft = fliplr(tempx((tempx > pmlin_minx-1e-8) & (tempx < pmlin_minx + hx+1e-8)));
xright = tempx((tempx > pmlin_maxx-hx-1e-8) & (tempx < pmlin_maxx +1e-8));

Nxl = length(xleft);
Nxr = length(xright);

wleft = hann(2*(Nxl-1)+1);
wleft = wleft(Nxl:end);

wright = hann(2*(Nxr-1)+1);
wright = wright(Nxr:end);

wind = ones(size(x(surf_node)));

for i = 1:Nxl
    ind = find((tempx > xleft(i)-1e-6) & (tempx < xleft(i)+1e-6));
    wind(ind) = wind(ind) * wleft(i);
end

for i = 1:Nxr
    ind = find((tempx > xright(i)-1e-6) & (tempx < xright(i)+1e-6));
    wind(ind) = wind(ind) * wright(i);
end

indl = find(tempx < pmlin_minx);
indr = find(tempx > pmlin_maxx);
wind(indl) = 0;
wind(indr) = 0;


% *************************************************************************
% Plot results and calculate element quality
% *************************************************************************
plot_flag = 1; % flag showing whether the mesh will be plotted or not
if plot_flag
    
    figure; hold on
    triplot(conn, x, y)
    axis equal tight;
    set(gcf,'Color',[1 1 1]);
    xlabel('x (m)');  ylabel('y (m)');
    
    hold on
    plot(x(pmlbin_node), y(pmlbin_node),'m.', 'LineWidth',1)
    plot(x(pmlbout_node), y(pmlbout_node),'k.', 'LineWidth',1)
    plot(x(pml_node), y(pml_node),'y.', 'LineWidth',1)
    
    plot(x(scatb_node), y(scatb_node),'r.', 'LineWidth',1)
    plot(xmid(scatb_elm), ymid(scatb_elm),'g.', 'LineWidth',1)
    plot(x(huygb_node), y(huygb_node),'k*', 'LineWidth',1)
    plot(xmid(huygb_elm), ymid(huygb_elm),'m*', 'LineWidth',1)
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
