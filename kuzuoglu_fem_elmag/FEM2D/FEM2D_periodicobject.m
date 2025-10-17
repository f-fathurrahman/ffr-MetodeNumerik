%% Name: FEM2D_periodic_object.m 
% 2D Scattering from an Array of Periodic Objects (Sec. 5.6.5)
%
% It solves the 2D Helmholtz equation with linear triangular elements 
% using periodic FEM formulation in a unit cell.
% Both TMz and TEz polarizations are considered.
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

global e0;
global nu0;
global mu0;
global k0;
global omg;

% Constants
%**********
c0  = 3*1e8;           % m/sec, velocity of light in free space
nu0 = 120*pi;          % ohm, intrinsic impedance of the free space
e0  = (1e-9)/(36*pi);  % F/m, permittivity of free space
mu0 = 4*pi*1e-7;       % H/m, permeability of free space

% Input parameters
%*****************
freq    = 300;          % MHz, frequency
freq    = freq*1e6;     % Hz, frequency
lambda0 = c0/freq;      % meter, wavelength
k0      = 2*pi/lambda0; % 1/meter, wavenumber
omg     = 2*pi*freq;    % rad/sec, radial frequency 

% Choose object type
scat_type = 'pec'; % 'pec' or 'diel' (PEC or dielectric object)

polrz = 'tm'; % 'TM' or 'TE' polarization

er = 4;  % dielectric constant of dielectric object (used if scat_type='diel')

phii = 30;           % degree, angle of incident field
phii = phii*pi/180;  % radian, angle of incident field

% Mesh Generation
%****************
delh = lambda0/40;   % element size

cellx = 1.0*lambda0;  % length of the cell along x
celly = 1.0*lambda0;  % length of the cell along y

scat_shape = 1; 
% 1: circular object 
% 2: half-circular object 
% 3: elliptical object 
% 4: polygonal object 
% 5: cross-section of conesphere object 
% 6: ogive object 

disp('The mesh is being created...');
if (scat_shape <= 2)
    rc = 0.3*lambda0;      % radius of circular or half-circular scatterer
    obj = struct('rc', rc,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape, 'cellx', cellx, 'celly', celly);  
 
elseif scat_shape == 3
    rex = 0.3*lambda0;      % x-dim of ellipse
    rey = 0.2*lambda0;      % y-dim of ellipse    
    obj = struct('rex', rex, 'rey', rey,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape, 'cellx', cellx, 'celly', celly);  
   
elseif scat_shape == 4
    svx = [0.3 0.2 0.2 -0.2 -0.2 -0.3 -0.3 0.3];   % Scatterer Vertex Coordinates (x)
    svy = [0.2 0.2 -0.2 -0.2 0.2 0.2 -0.3 -0.3]; % Scatterer Vertex Coordinates (y)
    obj = struct('svx', svx, 'svy', svy,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape, 'cellx', cellx, 'celly', celly);    
   % Vertices must be in the counter-clockwise direction!
        
elseif scat_shape == 5
    rc = 0.3*lambda0;     % radius of circular base
    hang = 45 *pi/180;    % half cone angle
    hh = rc/tan(hang);    % height of the cone
    obj = struct('rc', rc, 'hang', hang, 'hh', hh,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape, 'cellx', cellx, 'celly', celly);    

elseif scat_shape == 6
    hang = 37.5 *pi/180;          % angle
    hh = 0.3*lambda0;             % half height
    rc = hh / (1/tan(hang*0.5));  % half height
    obj = struct('rc', rc, 'hang', hang, 'hh', hh,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape, 'cellx', cellx, 'celly', celly);    
end

[conn,x,y,pmlbint_node,pmlboutt_node,pmlt_node,...
    pmlbinb_node,pmlboutb_node,pmlb_node,scatb_node,...
    bright_node,bleft_node,scatb_elm,scatin_node,scatin_elm] = meshgen_periodicobject(obj,freq);

disp('The mesh has been created!');

M = size(conn,1); % number of elements
N = length(x);    % number of nodes

% Material and source parameters
%*******************************
epsr = ones(M,1);    % relative permittivity of elements
epsr(scatin_elm) = er;  
mur = 1;             % relative permeability


%% ***********************************************************************
% MAIN BODY
%*************************************************************************

% Implement Locally-Conformal PML (LC-PML)
%*****************************************
pmlbint_x = x(pmlbint_node);    pmlbint_y = y(pmlbint_node);
pmlboutt_x = x(pmlboutt_node);  pmlboutt_y = y(pmlboutt_node);
pmlbinb_x = x(pmlbinb_node);    pmlbinb_y = y(pmlbinb_node);
pmlboutb_x = x(pmlboutb_node);  pmlboutb_y = y(pmlboutb_node);

disp('The LC-PML is being implemented...');
for i = 1:length(pmlt_node)
    in = pmlt_node(i);
    [x(in), y(in)] = lcpml(x(in), y(in), k0, pmlbint_x, pmlbint_y, pmlboutt_x, pmlboutt_y);
end    
for i = 1:length(pmlb_node)
    in = pmlb_node(i);
    [x(in), y(in)] = lcpml(x(in), y(in), k0, pmlbinb_x, pmlbinb_y, pmlboutb_x, pmlboutb_y);
end     
disp('The LC-PML has been implemented!');

% Matrix formation
%*****************
Ne = 3;           % number of nodes in each triangular element
A = sparse(N,N);  % initialize the global matrix to sparse zero matrix (size: NxN)
b = sparse(N,1);  % initialize the global b vector to sparse zero vector (size: Nx1)

uinc = exp(1j*k0*(real(x)*cos(phii)+real(y)*sin(phii))); % incident field at nodes

disp('The matrix is being constructed...');
% Start the matrix assembly process:
for e = 1:M  % loop over all elements
    if strcmpi(polrz, 'tm') % TM
        pe = 1/mur;
        qe = -k0^2*epsr(e);
    else % TE
        pe = 1/epsr(e);
        qe = -k0^2*mur;
    end

    [Ae, ~, ~] = element_matrix_tri2(e, x, y, conn, pe, pe, qe, 0);    % element matrix and b-vector
    %[Ae, ~, ~] = element_matrix_tri(e, x, y, conn, pe, pe, qe, 0);    % element matrix and b-vector
    
    for i = 1:Ne            % loop over the local nodes of each element
        ig = conn(e,i);     % global node corresponding to i
        for j = 1:Ne        % loop over the local nodes of each element
            jg = conn(e,j); % global node corresponding to j
            A(ig,jg) = A(ig,jg) + Ae(i,j); % fill the global matrix
        end
    end
end

if strcmpi(scat_type, 'diel')  % dielectric
    Au = A*uinc;
    b(scatin_node) = -1.0*Au(scatin_node);
end

disp('The matrix has been constructed!');

% Imposition of BCs
%**********************
if strcmpi(scat_type, 'pec') % PEC object
    disp('The BCs are being imposed...');

    if strcmpi(polrz, 'tm') % TMz
        BCnodes = [scatb_node; pmlboutt_node; pmlboutb_node];  % boundary node IDs
        BCvalues = [-uinc(scatb_node); zeros(size(pmlboutt_node)); zeros(size(pmlboutb_node))];  % function values at the boundary nodes        
        % the following is a fast implementation of Dirichlet BCs
        nodes = 1:N;
        nodes(BCnodes) = 0;
        othernodes = find(nodes ~= 0);
        Atemp = speye(N, N);
        Atemp(othernodes,:) = 0;
        A(BCnodes,:) = 0;
        A = A+Atemp;
        b(BCnodes) = BCvalues;
        
    else % TEz
        snodes = zeros(1,N);
        snodes(scatb_node) = 1;
        for i = 1:length(scatb_elm)
            e = scatb_elm(i);
            n1 = conn(e,1);  n2 = conn(e,2);  n3 = conn(e,3);
            
            if (snodes(n1) && snodes(n2) && ~snodes(n3))
                node1 = n1;  node2 = n2;  nodeout = n3;
            elseif (~snodes(n1) && snodes(n2) && snodes(n3))
                node1 = n2;  node2 = n3;  nodeout = n1;
            elseif (snodes(n1) && ~snodes(n2) && snodes(n3))
                node1 = n3;  node2 = n1;  nodeout = n2;
            else
                disp('Error in the BC implementation.');
                return;
            end
            
            [anx, any] = unit_vector(x, y, node1, node2, nodeout);
            anx = -anx; any = -any;
            
            delc = sqrt((x(node1)-x(node2))^2+(y(node1)-y(node2))^2);
            ym = 0.5*(y(node1)+y(node2));
            xm = 0.5*(x(node1)+x(node2));
            
            um = exp(1j*k0*(xm*cos(phii)+ym*sin(phii)));
            g = -1j*k0*um*(cos(phii)*anx + sin(phii)*any);
            term = g*delc*0.5;
            
            b(node1) = b(node1)+term;
            b(node2) = b(node2)+term;
        end
    end    
    disp('The BCs have been imposed!');
end

% Imposition of Periodic BCs
%*******************************
Nb = length(bright_node);
dx = abs(max(real(x)) - min(real(x)));

for ii = 1:Nb
    iright = bright_node(ii);
    
    ind = find(real(y(bleft_node)) == real(y(iright)));    
    ileft = bleft_node(ind);
    
    A(ileft, :) = A(ileft, :) + A(iright, :)*exp(-1j*k0*dx*cos(phii));
    b(ileft) = b(ileft) + b(iright)*exp(-1j*k0*dx*cos(phii));    
    
    A(iright, :) = 0;
    A(iright, iright) = 1;
    A(iright, ileft) = -exp(1j*k0*dx*cos(phii));
    b(iright) = 0;    
end

% Solution of the global matrix equation
%***************************************
disp('The matrix equation is being solved...');

u = A\b;         % sattered field values at the nodes
utot = u + uinc; % total field

disp('The matrix equation has been solved!');


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************

% Plot the scattered field
%*************************
[up, xp, yp] = create2darray(real(x), real(y), abs(u), delh, delh);

figure, clf, whitebg('white')
imagesc(xp, yp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Scattered Field (mag)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(x(scatb_node), y(scatb_node),'k.','linewidth',2)
plot(x(pmlbinb_node), y(pmlbinb_node),'k.','linewidth',1)
plot(x(pmlbint_node), y(pmlbint_node),'k.','linewidth',1)

% Plot the total field
%*********************
[up, xp, yp] = create2darray(real(x), real(y), abs(utot), delh, delh);

figure, clf, whitebg('white')
imagesc(xp, yp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Total Field (mag)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(x(scatb_node), y(scatb_node),'k.','linewidth',2)
plot(x(pmlbinb_node), y(pmlbinb_node),'k.','linewidth',1)
plot(x(pmlbint_node), y(pmlbint_node),'k.','linewidth',1)

disp('THE END !');

toc
