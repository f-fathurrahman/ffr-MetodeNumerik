%% Name: FEM2D_wg2.m 
% Scattering from an object within a parallel-plate waveguide (Sec. 5.6.8)
%
% It solves the 2D Helmholtz equation with linear triangular elements.
% Both TMz and TEz polarizations are considered.
% Reflection and transmission coefficients are determined. 
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

% Constants
%**********
c0  = 3*1e8;           % m/sec, velocity of light in free space
nu0 = 120*pi;          % ohm, intrinsic impedance of the free space
e0  = (1e-9)/(36*pi);  % F/m, permittivity of free space
mu0 = 4*pi*1e-7;       % H/m, permeability of free space

% Input parameters
%*****************
freq    = 2000;         % MHz, frequency
freq    = freq*1e6;     % Hz, frequency
lambda0 = c0/freq;      % meter, wavelength
k0      = 2*pi/lambda0; % 1/meter, wavenumber
omg     = 2*pi*freq;    % rad/sec, radial frequency 

polrz = 'tm'; % 'TM' or 'TE' polarization

er = 4;  % dielectric constant of dielectric object (used if scat_type='diel')
mur = 1; % relative permittivity of object (used if scat_type='diel')

% Mesh Generation
%****************
delh = 0.1/15;     % m, element size (lambda0/40)
dz = 0.5;          % m, length of the waveguide along z
dx = 0.1;          % m, length of the waveguide along x
scat_type = 'pec'; % 'pec' or 'diel' (PEC or dielectric object)
scat_shape = 2;    % 1: circular, 2: rectangular
rc = 0.03;         % m, radius of circular object
lz = 0.01/2;       % m, half-length of rectangular object along z
lx = 0.1/8;        % m, half-length of rectangular object along x
cx = dx/2;%lx/2;         % m, center of the object along x
pmldim = 0.05;     % m, distance btw inner and outer PML boundaries
dz1 = 0.05;        % m, z1 place for reflection coefficient measurement
dz2 = 0.45;        % m, z2 place for transmission coefficient measurement

obj = struct('dz', dz, 'dx', dx, 'scat_type', scat_type, ...
    'scat_shape', scat_shape, 'rc', rc, 'lz', lz, 'lx', lx, 'cx', cx,...
     'pmldim', pmldim,'dz1', dz1,'dz2', dz2, 'delh', delh);

disp('The mesh is being created...');
[conn,z,x,pmlbinr_node,pmlboutr_node,pmlr_node, ...
    pmlbinl_node,pmlboutl_node,pmll_node,scatb_node,btop_node, ...
    bbottom_node,scatin_node,scatin_elm,scatb_elm,dz1_node,dz2_node] ...
    = meshgen_wg2(obj);  
disp('The mesh has been created!');

M = size(conn,1); % number of elements
N = length(z);    % number of nodes

% Cutoff frequency and wavenumbers
%*********************************
nmode = 1;
kc = nmode*pi/dx;
fc = kc/(2*pi*sqrt(mu0*e0));
beta = sqrt(k0^2 - kc^2);

disp(['Cutoff freq (MHz) = ' sprintf('%.2f', fc*1e-6)]);
disp(['kc = ' sprintf('%.4f', kc)]);
disp(['k = ' sprintf('%.4f', k0)]);
disp(['beta = ' sprintf('%.4f', beta)]);

if imag(beta) ~= 0
    disp('Error: Beta is not real! Increase frequency!');
    return;
end

% Material and source parameters
%*******************************
epsr = ones(M,1);    % relative permittivity of elements
epsr(scatin_elm) = er;  


%% ***********************************************************************
% MAIN BODY
%*************************************************************************

% Implement Locally-Conformal PML (LC-PML)
%*****************************************
pmlbinr_z = z(pmlbinr_node);    pmlbinr_x = x(pmlbinr_node);
pmlboutr_z = z(pmlboutr_node);  pmlboutr_x = x(pmlboutr_node);
pmlbinl_z = z(pmlbinl_node);    pmlbinl_x = x(pmlbinl_node);
pmlboutl_z = z(pmlboutl_node);  pmlboutl_x = x(pmlboutl_node);

disp('The LC-PML is being implemented...');
for i = 1:length(pmlr_node)
    in = pmlr_node(i);
    [z(in), x(in)] = lcpml(z(in), x(in), k0, pmlbinr_z, pmlbinr_x, pmlboutr_z, pmlboutr_x);
end   
for i = 1:length(pmll_node)
    in = pmll_node(i);
    [z(in), x(in)] = lcpml(z(in), x(in), k0, pmlbinl_z, pmlbinl_x, pmlboutl_z, pmlboutl_x);
end
disp('The LC-PML has been implemented!');

% Matrix formation
%*****************
Ne = 3;           % number of nodes in each triangular element
A = sparse(N,N);  % initialize the global matrix to sparse zero matrix (size: NxN)
b = sparse(N,1);  % initialize the global b vector to sparse zero vector (size: Nx1)

if strcmpi(polrz, 'tm') % TMz
   uinc = sin(nmode*pi*real(x)/dx).*exp(-1j*beta*(real(z)-1j*abs(imag(z)))); % incident field at nodes
else
   uinc = cos(nmode*pi*real(x)/dx).*exp(-1j*beta*(real(z)-1j*abs(imag(z)))); % incident field at nodes
end

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

    [Ae, ~, ~] = element_matrix_tri2(e, z, x, conn, pe, pe, qe, 0);    % element matrix and b-vector
    %[Ae, ~, ~] = element_matrix_tri(e, z, x, conn, pe, pe, qe, 0);    % element matrix and b-vector
    
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

% Imposition of the BCs 
%**********************
disp('The BCs are being imposed...');

if strcmpi(polrz, 'tm') % TMz
    if strcmpi(scat_type, 'pec') % PEC object
        BCnodes = [scatb_node; btop_node; bbottom_node; pmlboutl_node; pmlboutr_node];  % boundary node IDs
        BCvalues = [-uinc(scatb_node); zeros(size(btop_node)); zeros(size(bbottom_node)); zeros(size(pmlboutl_node)); zeros(size(pmlboutr_node))];  % function values at the boundary nodes        
    else
        BCnodes = [btop_node; bbottom_node];  % boundary node IDs
        BCvalues = [zeros(size(btop_node)); zeros(size(bbottom_node))];  % function values at the boundary nodes
    end
    
    % the following is a fast implementation of the Dirichlet BCs
    nodes = 1:N;
    nodes(BCnodes) = 0;
    othernodes = find(nodes ~= 0);
    Atemp = speye(N, N);
    Atemp(othernodes,:) = 0;
    A(BCnodes,:) = 0;
    A = A+Atemp;
    b(BCnodes) = BCvalues;
end

if strcmpi(polrz, 'te') && strcmpi(scat_type, 'pec')  % TEz
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
        
        [anz, anx] = unit_vector(z, x, node1, node2, nodeout);
        anz = -anz; anx = -anx;
        
        delc = sqrt((z(node1)-z(node2))^2+(x(node1)-x(node2))^2);
        xm = 0.5*(x(node1)+x(node2));
        zm = 0.5*(z(node1)+z(node2));
        
        um = exp(-1j*beta*zm);
        g = -((-1j*beta*um*cos(nmode*pi*xm/dx))*anz + (-(nmode*pi/dx)*sin(nmode*pi*xm/dx)*um)*anx);
        term = g*delc*0.5;
        
        b(node1) = b(node1)+term;
        b(node2) = b(node2)+term;
    end
end
disp('The BCs have been imposed!');

% Solution of the global matrix equation
%***************************************
disp('The matrix equation is being solved...');

u = A\b;         % sattered field values at the nodes
utot = u + uinc; % total field

disp('The matrix equation has been solved!');


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************
R = max(abs(u(dz1_node)));    % reflection coefficient
T = max(abs(utot(dz2_node))); % transmission coefficient
disp(['R = ' sprintf('%.4f',full(R))]);
disp(['T = ' sprintf('%.4f',full(T))]);
disp(['Accuracy check: |R|^2+|T|^2 = ' sprintf('%.4f',R^2+T^2) ' =~ 1'])

% Plot the fields at dz1 and dz2
%*******************************
figure, clf, whitebg('white'), hold on
set(gcf, 'Color', [1 1 1])
plot(abs(u(dz1_node)),real(x(dz1_node)),'k', 'Linewidth', [2]); axis tight;
title('Reflected field @ dz1');
ylabel('x (m)');
xlabel('Field (mag)');

figure, clf, whitebg('white'), hold on
set(gcf, 'Color', [1 1 1])
plot(abs(utot(dz2_node)), real(x(dz2_node)),'k', 'Linewidth', [2]); axis tight;
title('Transmitted field @ dz2');
ylabel('x (m)');
xlabel('Field (mag)');

% Plot the scattered field
%*************************
[up, zp, xp] = create2darray(real(z), real(x), abs(u), delh, delh);
figure, clf, whitebg('white')
imagesc(zp, xp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('Scattered Field (mag)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(z(scatb_node), x(scatb_node),'k.','linewidth',3)

[up, zp, xp] = create2darray(real(z), real(x), real(u), delh, delh);
figure, clf, whitebg('white')
imagesc(zp, xp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('Scattered Field (real)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(z(scatb_node), x(scatb_node),'k.','linewidth',3)

% Plot the total field
%****************************************************
[up, zp, xp] = create2darray(real(z), real(x), abs(utot), delh, delh);
figure, clf, whitebg('white')
imagesc(zp, xp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('Total Field (mag)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(z(scatb_node), x(scatb_node),'k.','linewidth',3)

[up, zp, xp] = create2darray(real(z), real(x), real(utot), delh, delh);
figure, clf, whitebg('white')
imagesc(zp, xp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('Total Field (mag)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(z(scatb_node), x(scatb_node),'k.','linewidth',3)

disp('THE END !');
toc
