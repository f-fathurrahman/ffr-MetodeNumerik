%% Name: FEM2D_cap.m 
% 2D Parallel-Plate Capacitor Problem with Truncated Computational Domain (Sec. 5.5.1)
%
% It solves the 2D electrostatic problem with linear triangular or 
% bilinear quadrilateral elements.
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

%% ***********************************************************************
% PRE-PROCESSING
%*************************************************************************
clear all; close all;
tic

% Choose element type
%********************
elm_type = 1; % 1: triangular, 2: quadrilateral

% Mesh Generation
%****************
Lx = .05; % m, length of the computational domain along x axis
Ly = .05; % m, length of the computational domain along y axis

dx = .02; % m, length of the capacitor plate along x axis
dy = .01; % m, length of the dielectric region inside capacitor along y axis

delx = 0.0005; % m, increment along x axis
dely = delx;   % m, increment along y axis

% call mesh generation function
if elm_type == 1  % triangular
   [conn,x,y,xmid,ymid,outb_node,ptop_node,pbottom_node,d_elm] = meshgen_cap(Lx,Ly,dx,dy,delx,dely);
elseif elm_type == 2  % quadrilateral
   [conn,x,y,xmid,ymid,outb_node,ptop_node,pbottom_node,d_elm] = meshgenq_cap(Lx,Ly,dx,dy,delx,dely);        
else
   disp('Element type can only be triangular or quadriateral.');
   return;
end

M = size(conn,1); % number of elements
N = length(x);    % number of nodes

% Input parameters and constants
%*******************************
e0 = 8.85e-12; % free-space permittivity (F/m)
er = 1;        % dielectric constant inside the capacitor

% Boundary conditions
%********************
Vtop = 1;    % Dirichlet BC at the top plate
Vbottom = 0; % Dirichlet BC at the bottom plate

% Material and source parameters
%*******************************
eps = ones(M,1)*e0;          
eps(d_elm) = er*e0;  % permittivity (F/m) of elements


%% ***********************************************************************
% MAIN BODY
%*************************************************************************
if elm_type == 1  % triangular
   Ne = 3; % number of nodes in each triangular element
elseif elm_type == 2  % quadrilateral
   Ne = 4; % number of nodes in each quadrilateral element 
end

A = sparse(N,N);   % initialize the global matrix to sparse zero matrix (size: NxN)
b = sparse(N,1);   % initialize the global b vector to sparse zero vector (size: Nx1)
area = zeros(M,1); % initialize the area of elements(size: Mx1)

% Start the matrix assembly process:
%***********************************
for e = 1:M  % loop over all elements
    if elm_type == 1  % triangular
       [Ae, be, area(e)] = element_matrix_tri2(e, x, y, conn, eps(e), eps(e), 0, 0);   % element matrix and b-vector
%       [Ae, be, area(e)] = element_matrix_tri(e, x, y, conn, eps(e), eps(e), 0, 0);   % element matrix and b-vector       
    elseif elm_type == 2  % quadrilateral
        [Ae, be, area(e)] = element_matrix_quad2(e, x, y, conn, eps(e), eps(e), 0, 0); % element matrix and b-vector
%        [Ae, be, area(e)] = element_matrix_quad(e, x, y, conn, eps(e), eps(e), 0, 0); % element matrix and b-vector        
    end
    
    for i = 1:Ne            % loop over the local nodes of each element
        ig = conn(e,i);     % global node corresponding to i
        for j = 1:Ne        % loop over the local nodes of each element
            jg = conn(e,j); % global node corresponding to j
            A(ig,jg) = A(ig,jg) + Ae(i,j); % fill the global matrix
        end
        b(ig) = b(ig) + be(i); % fill the global b vector
    end
end

% Imposition of Dirichlet boundary conditions (first approach)
%*****************************************************************
BCnodes = [outb_node; ptop_node; pbottom_node];  % boundary node IDs
BCvalues = [zeros(size(outb_node)); Vtop*ones(size(ptop_node)); ...
            Vbottom*ones(size(pbottom_node))];  % function values at the boundary nodes

% the following is a fast implementation of Dirichlet BCs
nodes = 1:N;
nodes(BCnodes) = 0;
othernodes = find(nodes ~= 0);
Atemp = speye(N, N);
Atemp(othernodes,:) = 0;
A(BCnodes,:) = 0;
A = A+Atemp;
b(BCnodes) = BCvalues;

% Solution of the global matrix equation
%***************************************
V = A\b;     % potential values at the nodes


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************

% Electric field within each element
%***********************************
Ex = zeros(M,1); Ey = zeros(M,1);
if elm_type == 1  % triangular
    for e = 1:M  
        [Ex(e), Ey(e)] = Efield(e, x, y, conn, V);
    end
elseif elm_type == 2  % quadrilateral
    for e = 1:M  
        [Ex(e), Ey(e)] = Efieldq(e, x, y, conn, V);
    end
end
Emag = sqrt(Ex.^2 + Ey .^2); % magnitude of the E-field

% Capacitance
%************
%C = (1/(Vtop-Vbottom)^2)*sum(eps(d_elm).*(Ex(d_elm).^2 + Ey(d_elm).^2).*area(d_elm)); % F/m
C = (1/(Vtop-Vbottom)^2)*sum(eps.*(Ex.^2 + Ey.^2).*area); % F/m

disp(['C = ' num2str(C*1e12) ' pF/m'])

% Plot the potential
%*******************
[Vp, xp, yp] = create2darray(x, y, V, delx, dely);

figure, clf, whitebg('white')
imagesc(xp, yp, Vp);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Potential (V)');
set(gcf, 'Color', [1 1 1]);
hold on
% quiver(xmid, ymid, Ex, Ey, 'color', 'k'); % plot Efield vectors
plot(x(ptop_node), y(ptop_node),'k','linewidth',5)
plot(x(pbottom_node), y(pbottom_node),'k','linewidth',5)

figure, clf, whitebg('white')
contour(xp, yp, Vp, 'LevelStep',0.2)
colormap(jet)
axis equal tight;
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Potential (V)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(x(ptop_node), y(ptop_node),'k','linewidth',5)
plot(x(pbottom_node), y(pbottom_node),'k','linewidth',5)

% Plot the Efield
%****************
[Emagp, xp, yp] = create2darray(xmid, ymid, Emag, delx, dely);

figure, clf, whitebg('white')
imagesc(xp, yp, Emagp);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Electric Field Intensity (mag) (V/m)');
set(gcf, 'Color', [1 1 1]);
hold on
quiver(xmid, ymid, Ex, Ey, 'color', 'k'); % plot Efield vectors
plot(x(ptop_node), y(ptop_node),'k','linewidth',5)
plot(x(pbottom_node), y(pbottom_node),'k','linewidth',5)

toc
