function [conn,x,y,pmlbin_node,pmlbout_node,pml_node,scatb_node,...
    bright_node,bleft_node,scatb_elm] = meshgen_periodicgrating(obj,freq)
% Mesh generation function for "Scattering from a PEC Periodic Grating"
% using triangular elements
% Used in FEM2D_periodicgrating.m
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
% scatb_node : Array containing the nodes on the surface of the grating
% scatb_elm : Array containing the elements on the surface of the grating
% bright_node : Array containing the nodes on the right periodic boundary
% bleft_node : Array containing the nodes on the left periodic boundary
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

c0 = 3*1e8;        % m/sec, velocity of light in free space
lambda = c0/freq;  % meter, wavelength

delh = obj.delh;   % element size

svx = obj.svx;    % Scatterer Vertex Coordinates (x)
svy = obj.svy;    % Scatterer Vertex Coordinates (y)
svx = round(svx/delh)*delh; 
svy = round(svy/delh)*delh; 

fsdim = 1.0*lambda;      % distance btw scatterer and inner PML boundary
pmldim = 0.5*lambda;     % distance btw inner and outer PML boundaries
% fit to delh incrementation
fsdim = round(fsdim/delh)*delh; 
pmldim = round(pmldim/delh)*delh; 

pmlin_maxy = max(svy)+fsdim;      % inner PML boundary (max_y)
pmlout_maxy = pmlin_maxy+pmldim;  % outer PML boundary (max_y)

% *************************************************************************
% Create coordinates 
% *************************************************************************
Nx = round((max(svx)-min(svx))/delh)+1; 
Ny = round((pmlout_maxy-min(svy))/delh)+1; 

xt = linspace(min(svx), max(svx), Nx);
yt = linspace(min(svy), pmlout_maxy, Ny);

[xa, ya] = meshgrid(xt, yt);
xat = xa'; yat = ya';
x = xat(:); y = yat(:);
clear xat; clear yat;
    
svdist = sqrt((svx(2:end)-svx(1:end-1)).^2+(svy(2:end)-svy(1:end-1)).^2);
nd = round(svdist/delh)+1;
NV = length(svdist);

scatb_x = []; scatb_y = [];
for i = 1:NV
    tempx = linspace(svx(i), svx(i+1), nd(i));
    tempy = linspace(svy(i), svy(i+1), nd(i));
    scatb_x = [scatb_x; tempx'];
    scatb_y = [scatb_y; tempy'];
end
   
% *************************************************************************
threshold = delh*0.8;
[discard,d] = knnsearch([x y],[scatb_x scatb_y],'k',10);
ind = find(d < threshold);
discard = unique(discard(ind));
nodes0 = ones(1,length(x));
nodes0(discard) = 0;
nondiscard = find(nodes0 == 1);
x = x(nondiscard);
y = y(nondiscard);
x = cat(1, x, scatb_x);
y = cat(1, y, scatb_y);
  
coords = [x y];
coords = unique(coords,'rows');
x = coords(:,1);
y = coords(:,2);

% *************************************************************************
% Create mesh and connectivity matrix
% *************************************************************************
conn = delaunay(x, y);    

N = length(x);
M = size(conn,1);

% *************************************************************************
% Find discarded elements and modify mesh if scattterer is PEC
% (make inside of scatterer empty)
% *************************************************************************
xmid = mean(x(conn),2);
ymid = mean(y(conn),2);

INPO = inpolygon(xmid, ymid, [svx(1) svx svx(end)], [svy(1)-1 svy svy(end)-1]);
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

x = x(nondiscard_node);
y = y(nondiscard_node);

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

N = length(x);
M = size(conn,1);
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

xmid = mean(x(conn),2);
ymid = mean(y(conn),2);

% *************************************************************************
% Find scatterer boundary nodes
% *************************************************************************
[scatb_node,~] = knnsearch([x y],[scatb_x scatb_y],'k',1);
scatb_node = sort(unique(scatb_node));

% *************************************************************************
% Find PML boundary nodes
% *************************************************************************
pmlbin_node = find((y > pmlin_maxy-1e-8) & (y < pmlin_maxy+1e-8));
pmlbout_node = find((y > pmlout_maxy-1e-8) & (y < pmlout_maxy+1e-8));

% *************************************************************************
% Find PML inner nodes
% *************************************************************************
pml_node = find((y > pmlin_maxy+1e-8) & (y < pmlout_maxy+1e-8));

% *************************************************************************
% Find periodic boundary nodes
% *************************************************************************
bright_node = find((x > max(x)-1e-8) & (x < max(x)+1e-8));
bleft_node = find((x > min(x)-1e-8) & (x < min(x)+1e-8));

% *************************************************************************
% Find scatterer boundary elements
% *************************************************************************
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
    
% *************************************************************************
% Plot results
% *************************************************************************
plot_flag = 1; % flag showing whether the mesh will be plotted or not
if plot_flag
    figure; hold on
    triplot(conn, x, y)
    axis equal tight;
    set(gcf,'Color',[1 1 1]);
    xlabel('x (m)');  ylabel('y (m)');
    
    plot(x(scatb_node), y(scatb_node),'r.', 'LineWidth',1.5)
    plot(x(pmlbout_node), y(pmlbout_node),'k.', 'LineWidth',1.5)
    plot(x(pmlbin_node), y(pmlbin_node),'m.', 'LineWidth',1.5)
    plot(x(pml_node), y(pml_node), 'y.', 'LineWidth',1.5)
    plot(x(bright_node), y(bright_node), 'go', 'LineWidth',1.5)
    plot(x(bleft_node), y(bleft_node), 'co', 'LineWidth',1.5)  
    plot(xmid(scatb_elm), ymid(scatb_elm), 'g.')
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
