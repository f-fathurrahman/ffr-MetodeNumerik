function [conn,x,y,wgb_node] = meshgen_wg4(obj)
% Mesh generation function for "Propagation inside a 3-D waveguide with
% uniform cross section" using triangular elements
% Used in FEM2D_wg4.m
%
% INPUT:
% obj  : object structure (see the content of the structure in the main file)
% OUTPUT:
% conn : connectivity matrix (size: Mx3)
% x,y  : x and y coordinates of all nodes in the mesh (each size: Nx1)
% wgb_node : Array containing the nodes on the surface of the waveguide
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
delh = obj.delh;         % element size
wg_shape = obj.wg_shape; % 1: circular, 2: rectangular
rc = obj.rc;             % radius of circular scatterer
lx = obj.lx;             % half-length of the rectangle along x
ly = obj.ly;             % half-length of the rectangle along y

% fit to delh incrementation
rc = round(rc/delh)*delh; 
lx = round(lx/delh)*delh; 
ly = round(ly/delh)*delh; 

% *************************************************************************
% Create coordinates 
% *************************************************************************
if wg_shape == 1
    x = 0; y = 0; % add center to mesh
    for r = delh:delh:rc
        % create a single circle with radius r
        N = floor(pi/(2*asin(delh*0.5/r)));
        alfa = cat(2, (pi/N)*(0:N), -fliplr((pi/N)*(1:N-1)));
        x0 = r*cos(alfa)'; % x-coordinates of the circle
        y0 = r*sin(alfa)'; % y-ccordinates of the circle
        
        x = cat(1, x, x0);
        y = cat(1, y, y0);
        
        if (r < rc+1e-8) && (r > rc-1e-8)
            wgb_x = x0;
            wgb_y = y0;
        end
    end
    
else
    Nx = round(2*lx/delh)+1;
    Ny = round(2*ly/delh)+1;
    
    xt = linspace(-lx, lx, Nx);
    yt = linspace(-ly, ly, Ny);
    
    [xa, ya] = meshgrid(xt, yt);
    xat = xa'; yat = ya';
    x = xat(:); y = yat(:);
    clear xat; clear yat;
end

% *************************************************************************
% Create mesh and connectivity matrix
% *************************************************************************
conn = delaunay(x, y);    

N = length(x);
M = size(conn,1);

disp(['Nelement = ' sprintf('%d', M)]);
disp(['Nnode = ' sprintf('%d', N)]);

% *************************************************************************
% Find boundary nodes
% *************************************************************************
if wg_shape == 1
    [wgb_node,~] = knnsearch([x y],[wgb_x wgb_y],'k',1);
    wgb_node = sort(unique(wgb_node));
    
else
    wgb_node = find(((x < lx+1e-8) & (x > lx-1e-8) & ...
        (y < ly+1e-8) & (y > -ly-1e-8)) | ...
        ((x < -lx+1e-8) & (x > -lx-1e-8) & ...
        (y < ly+1e-8) & (y > -ly-1e-8)) | ...
        ((y < ly+1e-8) & (y > ly-1e-8) & ...
        (x < lx+1e-8) & (x > -lx-1e-8)) | ...
        ((y < -ly+1e-8) & (y > -ly-1e-8) & ...
        (x < lx+1e-8) & (x > -lx-1e-8)));
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
    
    plot(x(wgb_node), y(wgb_node),'r.', 'LineWidth',1.5)    
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
