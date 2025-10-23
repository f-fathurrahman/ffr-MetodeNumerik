function [conn,x,y,pmlbint_node,pmlboutt_node,pmlt_node,...
    pmlbinb_node,pmlboutb_node,pmlb_node,scatb_node,...
    bright_node,bleft_node,scatb_elm,scatin_node,scatin_elm] = meshgen_periodicobject(obj,freq)
% Mesh generation function for "Scattering from an Array of Periodic
% Objects" using triangular elements
% Used in FEM2D_periodicobject.m
%
% INPUT:
% obj  : object structure (see the content of the structure in the main file)
% freq : frequency (Hz)
% OUTPUT:
% conn : connectivity matrix (size: Mx3)
% x,y  : x and y coordinates of all nodes in the mesh (each size: Nx1)
% pmlbint_node : Array containing the nodes on the inner PML boundary (top)
% pmlboutt_node : Array containing the nodes on the outer PML boundary (top)
% pmlt_node : Array containing the nodes inside the PML region (top)
% pmlbinb_node : Array containing the nodes on the inner PML boundary (bottom)
% pmlboutb_node : Array containing the nodes on the outer PML boundary (bottom)
% pmlb_node : Array containing the nodes inside the PML region (bottom)
% scatb_node : Array containing the nodes on the surface of the scatterer
% scatb_elm : Array containing the elements on the surface of the scatterer
% scatin_node : Array containing the nodes inside the scatterer (empty for PEC scatterer)
% scatin_elm : Array containing the elements inside the scatterer (empty for PEC scatterer)
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

delh = obj.delh;    % element size
cellx = obj.cellx;  % length of the cell along x
celly = obj.celly;  % length of the cell along y

pmldim = 0.5*lambda;     % distance btw inner and outer PML boundaries
pmldim = round(pmldim/delh)*delh;  % fit to delh incrementation

scat_type = obj.scat_type;   % 'pec' or 'diel' (PEC or dielectric object)
scat_shape = obj.scat_shape; % 1,2,3,4 or 5

if scat_shape == 1 % circular object with rectangular domain
    rc = obj.rc;      % radius of circular scatterer
    rc = round(rc/delh)*delh; % fit to delh incrementation
elseif scat_shape == 2 % half-circular object with rectangular domain
    rc = obj.rc;      % radius of half-circular scatterer
    rc = round(rc/delh)*delh; % fit to delh incrementation
elseif scat_shape == 3 % elliptical object with rectangular domain 
    rex = obj.rex;      % x-dim of ellipse
    rey = obj.rey;      % y-dim of ellipse    
    rex = round(rex/delh)*delh; % fit to delh incrementation
    rey = round(rey/delh)*delh; % fit to delh incrementation    
elseif scat_shape == 4 % polygonal object with rectangular domain
    svx = obj.svx;    % Scatterer Vertex Coordinates (x)
    svy = obj.svy;    % Scatterer Vertex Coordinates (y)
    svx = round(svx/delh)*delh; 
    svy = round(svy/delh)*delh; 
elseif scat_shape == 5 % cross-section of conesphere object with rectangular domain    
    rc = obj.rc;          % radius of circular base
    hang = obj.hang;      % half cone angle
    hh = obj.hh;          % height of the cone
    rc = round(rc/delh)*delh;     % fit to delh incrementation
    hang = round(hang/delh)*delh; % fit to delh incrementation
    hh = round(hh/delh)*delh;     % fit to delh incrementation  
elseif scat_shape == 6  % ogive object with rectangular domain   
    rc = obj.rc;          % radius 
    hang = obj.hang;      % half angle
    hh = obj.hh;          % half height
    rc = round(rc/delh)*delh;     % fit to delh incrementation
    hang = round(hang/delh)*delh; % fit to delh incrementation
    hh = round(hh/delh)*delh;     % fit to delh incrementation  
end

%
if scat_shape == 1
    N = floor(pi/(2*asin(delh*0.5/rc)));
    alfa = cat(2, (pi/N)*(0:N), -fliplr((pi/N)*(1:N-1)));
    scatb_x = rc*cos(alfa)'; % x-coordinates of the circle
    scatb_y = rc*sin(alfa)'; % y-ccordinates of the circle

elseif scat_shape == 2    
    N = floor(pi/(2*asin(delh*0.5/rc)));
    alfa = (pi/N)*(0:N);
    scatb_x = rc*cos(alfa)'; % x-coordinates of the circle
    scatb_y = rc*sin(alfa)'; % y-ccordinates of the circle    
    
    nc = round(2*rc/delh)+1;
    tempr = linspace(-rc, rc, nc);
    temprx = tempr(2:end-1);
    tempry = zeros(size(temprx));
    scatb_x = [scatb_x; temprx'];
    scatb_y = [scatb_y; tempry'];
    
elseif scat_shape == 3    
    if rex >= rey
        xtemp2 = linspace(rex, rex-delh, 10);
        ytemp2 = rey*sqrt(1-xtemp2.^2/rex^2);
        threshold3 = 0.5;
        ymid = (ytemp2(1)+ytemp2(end))*threshold3;
        ytemp1 = abs(ytemp2-ymid);
        [~,ind] = min(ytemp1);
        delsm = rex-xtemp2(ind);
        
        xtemp1 = linspace(rex, -rex, round(2*rex/delh)+1);
        if ind == 1
            xtemp = xtemp1;
        else
            xtemp = cat(2, xtemp1(1), rex-delsm, xtemp1(2:end-1), -rex+delsm, xtemp1(end));
        end
        ytemp = rey*sqrt(1-xtemp.^2/rex^2);
        tempx = cat(2, xtemp, fliplr(xtemp(2:end-1)));
        tempy = cat(2, ytemp, -fliplr(ytemp(2:end-1)));
        clear xtemp; clear ytemp;
                
    else  %rex < rey,
        ytemp2 = linspace(rey, rey-delh, 10);
        xtemp2 = rex*sqrt(1-ytemp2.^2/rey^2);
        threshold3 = 0.5;
        xmid = (xtemp2(1)+xtemp2(end))*threshold3;
        xtemp1 = abs(xtemp2-xmid);
        [~,ind] = min(xtemp1);
        delsm = rey-ytemp2(ind);
        
        ytemp1 = linspace(rey, -rey, round(2*rey/delh)+1);
        if ind == 1
            ytemp = ytemp1;
        else
            ytemp = cat(2, ytemp1(1), rey-delsm, ytemp1(2:end-1), -rey+delsm, ytemp1(end));
        end
        xtemp = rex*sqrt(1-ytemp.^2/rey^2);
        tempy = cat(2, ytemp, fliplr(ytemp(2:end-1)));
        tempx = cat(2, xtemp, -fliplr(xtemp(2:end-1)));
        clear xtemp; clear ytemp;
        
    end
    scatb_x = tempx';
    scatb_y = tempy';          

elseif scat_shape == 4
    svx = [svx svx(1)];
    svy = [svy svy(1)];
    
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
   
elseif scat_shape == 5
    N = floor(pi/(2*asin(delh*0.5/rc)));
    alfa = -fliplr((pi/N)*(1:N-1));
    scatb_x = rc*cos(alfa)'; % x-coordinates of the circle
    scatb_y = rc*sin(alfa)'; % y-ccordinates of the circle
    
    disti = sqrt((rc)^2 + (hh)^2);
    nd = round(disti/delh)+1;    
    tempx1 = linspace(rc, 0, nd)';
    tempy1 = linspace(0, hh, nd)';
    tempx2 = linspace(0, -rc, nd)';
    tempy2 = linspace(hh, 0, nd)';
    
    scatb_x = cat(1, tempx1, tempx2(2:end), scatb_x);
    scatb_y = cat(1, tempy1, tempy2(2:end), scatb_y);

elseif scat_shape == 6
    cosa = cos(hang);
    minuscosa = 1 - cosa;
    sina = sin(hang);
    
    nx = round(2*hh/delh)+1;
    scatb_x = linspace(hh, -hh, nx)';
    fx = sqrt(1-(scatb_x/hh).^2*sina^2)-cosa;
    scatb_y = rc*fx ./ minuscosa;
    
    scatb_x = cat(1, scatb_x, flipud(scatb_x(2:end-1)));
    scatb_y = cat(1, scatb_y, -1.0*scatb_y(2:end-1));    
end

errorflag = 0;
if ((max(scatb_x)-min(scatb_x)) >= cellx) || ((max(scatb_y)-min(scatb_y)) >= celly) 
    error('ERROR: Object dimensions cannot be larger than the cell dimensions!');
end

scatb_x = scatb_x + cellx/2;
scatb_y = scatb_y + celly/2;

pmlin_maxy = celly;  % inner PML boundary (max_y)
pmlin_miny = 0;      % inner PML boundary (min_y)
pmlout_maxy = pmlin_maxy+pmldim;  % outer PML boundary (max_y)
pmlout_miny = pmlin_miny-pmldim;  % outer PML boundary (min_y)

% *************************************************************************
% Create coordinates 
% *************************************************************************
Nx = round(cellx/delh)+1; 
Ny = round((pmlout_maxy-pmlout_miny)/delh)+1; 

xt = linspace(0, cellx, Nx);
yt = linspace(pmlout_miny, pmlout_maxy, Ny);

[xa, ya] = meshgrid(xt, yt);
xat = xa'; yat = ya';
x = xat(:); y = yat(:);
clear xat; clear yat;
      
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

if ((scat_shape ~= 4)) && strcmpi(scat_type, 'pec')    % PEC  
    [IN, ON] = inpolygon(x, y, scatb_x, scatb_y);
    discard_node = find((IN == 1) & (ON == 0));
    
    nodes0 = ones(1,length(x));
    nodes0(discard_node) = 0;
    nondiscard_node = find(nodes0 ~= 0);
    x = x(nondiscard_node);
    y = y(nondiscard_node);    
    x = [x;mean(scatb_x)];
    y = [y;mean(scatb_y)];
end

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
if strcmpi(scat_type, 'pec')    % PEC
    xmid = mean(x(conn),2);
    ymid = mean(y(conn),2);
    
    INPO = inpolygon(xmid, ymid, scatb_x, scatb_y);
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
end

N = length(x);
M = size(conn,1);
disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

xmid = mean(x(conn),2);
ymid = mean(y(conn),2);

% *************************************************************************
% Find nodes and elements inside the scatterer if dielectric
% *************************************************************************
scatin_elm = [];  scatin_node = [];
if strcmpi(scat_type, 'diel') % diel
    INPO = inpolygon(xmid,ymid, scatb_x, scatb_y);
    scatin_elm = find(INPO == 1);
    
    INPO = inpolygon(x,y, scatb_x, scatb_y);
    scatin_node = find(INPO == 1);
end

% *************************************************************************
% Find scatterer boundary nodes
% *************************************************************************
[scatb_node,~] = knnsearch([x y],[scatb_x scatb_y],'k',1);
scatb_node = sort(unique(scatb_node));

% *************************************************************************
% Find PML boundary nodes (top and bottom)
% *************************************************************************
pmlbint_node = find((y > pmlin_maxy-1e-8) & (y < pmlin_maxy+1e-8));
pmlboutt_node = find((y > pmlout_maxy-1e-8) & (y < pmlout_maxy+1e-8));

pmlbinb_node = find((y > pmlin_miny-1e-8) & (y < pmlin_miny+1e-8));
pmlboutb_node = find((y > pmlout_miny-1e-8) & (y < pmlout_miny+1e-8));

% *************************************************************************
% Find PML inner nodes
% *************************************************************************
pmlt_node = find((y > pmlin_maxy+1e-8) & (y < pmlout_maxy+1e-8));
pmlb_node = find((y > pmlout_miny-1e-8) & (y < pmlin_miny-1e-8));

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

x = x-mean(scatb_x);
y = y-mean(scatb_y);
xmid = mean(x(conn),2);
ymid = mean(y(conn),2);

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
    plot(x(pmlboutt_node), y(pmlboutt_node),'k.', 'LineWidth',1.5)
    plot(x(pmlbint_node), y(pmlbint_node),'m.', 'LineWidth',1.5)
    plot(x(pmlt_node), y(pmlt_node), 'y.', 'LineWidth',1.5)
    plot(x(pmlboutb_node), y(pmlboutb_node),'k.', 'LineWidth',1.5)
    plot(x(pmlbinb_node), y(pmlbinb_node),'m.', 'LineWidth',1.5)
    plot(x(pmlb_node), y(pmlb_node), 'y.', 'LineWidth',1.5)    
    plot(x(bright_node), y(bright_node), 'go', 'LineWidth',1.5)
    plot(x(bleft_node), y(bleft_node), 'co', 'LineWidth',1.5)  
    plot(xmid(scatb_elm), ymid(scatb_elm), 'g.')
    plot(x(scatin_node), y(scatin_node), 'g.', 'LineWidth',1.5)
    plot(xmid(scatin_elm), ymid(scatin_elm), 'y.')    
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
