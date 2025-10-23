function [conn,x,y,pmlbin_node,pmlbout_node,pml_node,...
      scatb_node,scatin_node,scatin_elm,scatb_elm,huygb_node,huygb_elm] ...
      = meshgen_scat_circle(obj,freq)
% Mesh generation function for "Scattering from a PEC or Dielectric Object"
% using triangular elements
% Used for scat_shape = 1 (circular object with circular computational domain)
% Used in FEM2D_scat.m
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
% scatb_node : Array containing the nodes on the surface of the scatterer
% scatin_node : Array containing the nodes inside the scatterer (empty for PEC scatterer)
% scatin_elm : Array containing the elements inside the scatterer (empty for PEC scatterer)
% scatb_elm : Array containing the elements on the surface of the scatterer
% huygb_node : Array containing the nodes on the surface of the Huygens boundary
% huygb_elm : Array containing the elements on the surface of the Huygens boundary
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

delh = obj.delh;             % element size
scat_type = obj.scat_type;   % 'pec' or 'diel' (PEC or dielectric object)
rc = obj.rc;                 % radius of circular scatterer

fsdim = 0.5*lambda;      % distance btw scatterer and inner PML boundary
pmldim = 0.5*lambda;     % distance btw inner and outer PML boundaries
huygdim = delh;          % distance btw scatterer and Huygens boundary (used if scat_type='diel')

rpmlin = rc+fsdim;       % radius of inner PML boundary
rpmlout = rpmlin+pmldim; % radius of outer PML boundary
rhuyg = rc+huygdim;      % radius of Huygens boundary

% fit to delh incrementation
rc = round(rc/delh)*delh; 
rpmlout = round(rpmlout/delh)*delh; 
rpmlin = round(rpmlin/delh)*delh; 
rhuyg = round(rhuyg/delh)*delh; 

% *************************************************************************
% Create coordinates layer by layer
% *************************************************************************
x = []; y = [];

if strcmpi(scat_type, 'diel') % dielectric
    x = 0; y = 0; % add center to mesh
end

for r = delh:delh:rpmlout
    if (strcmpi(scat_type, 'diel')) || ((strcmpi(scat_type, 'pec')) && (r >= rc))
        % create a single circle with radius r
        N = floor(pi/(2*asin(delh*0.5/r)));
        alfa = cat(2, (pi/N)*(0:N), -fliplr((pi/N)*(1:N-1)));
        x0 = r*cos(alfa); % x-coordinates of the circle
        y0 = r*sin(alfa); % y-ccordinates of the circle
        
        x = cat(2, x, x0);
        y = cat(2, y, y0);
        
        if (r < rc+1e-8) && (r > rc-1e-8)
            scatb_x = x0;
            scatb_y = y0;
        elseif (r < rhuyg+1e-8) && (r > rhuyg-1e-8)
            huygb_x = x0;
            huygb_y = y0;
        elseif (r < rpmlin+1e-8) && (r > rpmlin-1e-8)
            pmlbin_x = x0;
            pmlbin_y = y0;
        elseif (r < rpmlout+1e-8) && (r > rpmlout-1e-8)
            pmlbout_x = x0;
            pmlbout_y = y0;
        end
    end
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
    discard_elm = find(sqrt(xmid.^2 + ymid.^2) < rc-1e-6);    
    
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

x = x'; y = y';
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
[scatb_node,~] = knnsearch([x y],[scatb_x' scatb_y'],'k',1);
scatb_node = sort(unique(scatb_node));

% *************************************************************************
% Find Huygens boundary nodes
% *************************************************************************
[huygb_node,~] = knnsearch([x y],[huygb_x' huygb_y'],'k',1);
huygb_node = sort(unique(huygb_node));

% *************************************************************************
% Find PML boundary nodes
% *************************************************************************
[pmlbin_node,~] = knnsearch([x y],[pmlbin_x' pmlbin_y'],'k',1);
pmlbin_node = sort(unique(pmlbin_node));

[pmlbout_node,~] = knnsearch([x y],[pmlbout_x' pmlbout_y'],'k',1);
pmlbout_node = sort(unique(pmlbout_node));

% *************************************************************************
% Find PML inner nodes
% *************************************************************************
INPO = inpolygon(x, y, pmlbin_x, pmlbin_y);
temp_node = find(INPO == 1);      

nodes0 = ones(1,N); 
nodes0(temp_node) = 0;
pml_node = find(nodes0 == 1);

% *************************************************************************
% Find scatterer or huygens boundary elements
% *************************************************************************
scatb_elm = []; huygb_elm = [];
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
        
        if (count >= 2)
            scatb_elm(ii) = e;
            ii = ii+1;
        end
    end
    
else % dielectric
    nodesall = zeros(1, N);
    nodesall(huygb_node) = 1;
    
    ii = 1;
    for e = 1:M
        
        count = 0;
        for i = 1:3
            if (nodesall(conn(e,i)) == 1)
                count = count+1;
            end
        end
        
        if (count >= 2)
            huygb_elm(ii) = e;
            ii = ii+1;
        end
    end
    
    INPO = inpolygon(xmid, ymid, x(huygb_node), y(huygb_node));
    discard_elm = find(INPO == 1);
    elms0 = zeros(1,M);
    elms0(huygb_elm) = 1;
    elms0(discard_elm) = 0;
    huygb_elm = find(elms0 == 1);
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
    
    plot([x(scatb_node); x(scatb_node(1))] , [y(scatb_node); y(scatb_node(1))],'r', 'LineWidth',1.5)
    plot([x(pmlbout_node); x(pmlbout_node(1))] , [y(pmlbout_node); y(pmlbout_node(1))],'k', 'LineWidth',1.5)
    plot([x(pmlbin_node); x(pmlbin_node(1))] , [y(pmlbin_node); y(pmlbin_node(1))],'m', 'LineWidth',1.5)
    plot(x(pml_node), y(pml_node), 'y.', 'LineWidth',1.5)
    plot(x(scatin_node), y(scatin_node), 'g.', 'LineWidth',1.5)
    plot(xmid(scatin_elm), ymid(scatin_elm), 'y.')
    plot(xmid(scatb_elm), ymid(scatb_elm), 'g*')
    plot([x(huygb_node); x(huygb_node(1))] , [y(huygb_node); y(huygb_node(1))],'g', 'LineWidth',1.5)
    plot(xmid(huygb_elm), ymid(huygb_elm), 'c*')
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
