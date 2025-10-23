%% Name: FEM2D_wg1.m 
% Visualization of modes in a parallel-plate waveguide and 2-D horn (Sec. 5.6.7)
%
% It solves the 2D Helmholtz equation with linear triangular elements.
% Both TMz and TEz polarizations are considered.
%
% It includes 4 examples:
% Example 1: Parallel plate waveguide with step discontinuity
% Example 2: Parallel plate waveguide mounted on ground
% Example 3: Horn mounted on ground
% Example 4: Horn
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
nu0 = 120*pi;          % ohm, intrinsic impedance of free space
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
nmode = 1;    % mode

example = 1; % enter 1, 2, 3 or 4
% Example 1: Parallel plate waveguide with step discontinuity
% Example 2: Parallel plate waveguide mounted on ground
% Example 3: Horn mounted on ground
% Example 4: Horn

% Mesh Generation
%****************
if example == 1 % Parallel plate waveguide with step discontinuity
    delh = 0.1/40;  % m, element size (lambda0/40)
    dx = 0.1;       % m, length of the waveguide along x
    dz = 0.5;       % m, length of the waveguide along z
    discontinuity_exist = 1; % flag showing whether the discontinuity exists or not
    lx = 0.05;      % m, length of the discontinuity along x
    lz = 0.1;       % m, length of the discontinuity along z
    pmldim = 0.1/2; % m, distance btw inner and outer PML boundaries (lambda0/4)
    
    obj = struct('dz', dz, 'dx', dx, 'discontinuity_exist', discontinuity_exist, ...
        'lz', lz, 'lx', lx, 'pmldim', pmldim, 'delh', delh);
    
    if (lx >= dx) || (lz >= dz)
        disp('ERROR: The dimensions of the discontinuity cannot be larger than the dimensions of the waveguide!')
        return;
    end
    
    disp('The mesh is being created...');
    [conn,z,x,pmlbin_node,pmlbout_node,pml_node,sourceb_node, ...
        btop_node,bbottom_node] = meshgen_wg1a(obj);
    disp('The mesh has been created!');
    
elseif example == 2 % Parallel plate waveguide mounted on ground
    delh = 0.1/20;  % m, element size (lambda0/40)
    dx = 0.1;       % m, length of the waveguide along x
    dz = 0.1;       % m, length of the waveguide along z
    lx = 0.1;       % m, length of the ground along x
    pmldim = 0.1/2; % m, distance btw inner and outer PML boundaries (lambda0/4)
    fsdim = 0.1;    % m, length of the free space region along z (btw ground and inner PML boundary)
    
    obj = struct('dz', dz, 'dx', dx, 'lx', lx, ...
        'pmldim', pmldim, 'fsdim', fsdim, 'delh', delh);
    
    disp('The mesh is being created...');
    [conn,z,x,pmlbin_node,pmlbout_node,pml_node,sourceb_node, ...
        btop_node,bbottom_node] = meshgen_wg1b(obj);
    disp('The mesh has been created!');
    
elseif example == 3 % Horn mounted on ground
    delh = 0.1/20;  % m, element size (lambda0/40)
    dx = 0.1;       % m, length of the waveguide along x
    dz = 0.05;      % m, length of the waveguide along z
    lx = 0.1;       % m, length of the ground along x
    ax = 0.05;      % m, length of the tilted edge of the horn along x
    az = 0.1;       % m, length of the tilted edge of the horn along z
    pmldim = 0.1/2; % m, distance btw inner and outer PML boundaries (lambda0/4)
    fsdim = 0.1;    % m, length of the free space region along z (btw ground and inner PML boundary)
    
    obj = struct('dz', dz, 'dx', dx, 'lx', lx, 'az', az, 'ax', ax, ...
        'pmldim', pmldim, 'fsdim', fsdim, 'delh', delh);
    
    disp('The mesh is being created...');
    [conn,z,x,pmlbin_node,pmlbout_node,pml_node,sourceb_node, ...
        btop_node,bbottom_node] = meshgen_wg1c(obj);
    disp('The mesh has been created!');
    
elseif example == 4 % Horn
    delh = 0.1/20;  % m, element size (lambda0/40)
    dx = 0.1;       % m, length of the waveguide along x
    dz = 0.1;      % m, length of the waveguide along z
    ax = 0.05;      % m, length of the tilted edge of the horn along x
    az = 0.1;       % m, length of the tilted edge of the horn along z
    pmldim = 0.1/2; % m, distance btw inner and outer PML boundaries (lambda0/4)
    fsdim = 0.1;    % m, length of the free space region along z (btw ground and inner PML boundary)
    
    obj = struct('dz', dz, 'dx', dx, 'az', az, 'ax', ax, ...
        'pmldim', pmldim, 'fsdim', fsdim, 'delh', delh);
    
    disp('The mesh is being created...');
    [conn,z,x,pmlbin_node,pmlbout_node,pml_node,sourceb_node, ...
        btop_node,bbottom_node] = meshgen_wg1d(obj);
    disp('The mesh has been created!');
    
end

M = size(conn,1); % number of elements
N = length(z);    % number of nodes

% Cutoff frequency and wavenumbers
%*********************************
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

if example == 1
   kcd = nmode*pi/lx;
   fcd = kcd/(2*pi*sqrt(mu0*e0));
   disp(['Cutoff freq of the discontinuity (MHz) = ' sprintf('%.2f', fcd*1e-6)]);
end

% Material and source parameters
%*******************************
epsr = ones(M,1); % relative permittivity (dielectric constant) of elements
mur = 1;          % relative permeability


%% ***********************************************************************
% MAIN BODY
%*************************************************************************

% Implement Locally-Conformal PML (LC-PML)
%*****************************************
pmlbin_z = z(pmlbin_node);    pmlbin_x = x(pmlbin_node);
pmlbout_z = z(pmlbout_node);  pmlbout_x = x(pmlbout_node);

disp('The LC-PML is being implemented...');
for i = 1:length(pml_node)
    in = pml_node(i);
    [z(in), x(in)] = lcpml(z(in), x(in), k0, pmlbin_z, pmlbin_x, pmlbout_z, pmlbout_x);
end       
disp('The LC-PML has been implemented!');

% Matrix formation
%*****************
Ne = 3;           % number of nodes in each triangular element
A = sparse(N,N);  % initialize the global matrix to sparse zero matrix (size: NxN)
b = sparse(N,1);  % initialize the global b vector to sparse zero vector (size: Nx1)

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

disp('The matrix has been constructed!');

% Imposition of BCs
%**********************
disp('The BCs are being imposed...');

if strcmpi(polrz, 'tm') % TMz
    term = sin(nmode*pi*x(sourceb_node)/dx).*exp(-1j*beta*z(sourceb_node));
    BCnodes = [sourceb_node; btop_node; bbottom_node];  % boundary node IDs
    BCvalues = [term; zeros(size(btop_node)); zeros(size(bbottom_node))];  % function values at the boundary nodes
else % TEz
    term = cos(nmode*pi*x(sourceb_node)/dx).*exp(-1j*beta*z(sourceb_node));
    BCnodes = sourceb_node;  % boundary node IDs
    BCvalues = term;  % function values at boundary nodes
end
% the following is a fast implementation of Dirichlet BCs
nodes = 1:N;
nodes(BCnodes) = 0;
othernodes = find(nodes ~= 0);
Atemp = speye(N, N);
Atemp(othernodes,:) = 0;
A(BCnodes,:) = 0;
A = A+Atemp;
b(BCnodes) = BCvalues;

disp('The BCs have been imposed!');

% Solution of the global matrix equation
%***************************************
disp('The matrix equation is being solved...');

u = A\b;  % field values at the nodes

disp('The matrix equation has been solved!');


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************

% Plot the field
%***************
[up, zp, xp] = create2darray(real(z), real(x), abs(u), delh, delh);
figure, clf, whitebg('white')
imagesc(zp, xp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('Field (mag)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(real(z(btop_node)), real(x(btop_node)),'k.','linewidth',3)
plot(real(z(bbottom_node)), real(x(bbottom_node)),'k.','linewidth',3)

figure, clf, whitebg('white')
contour(zp, xp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('Field (mag)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(real(z(btop_node)), real(x(btop_node)),'k.','linewidth',3)
plot(real(z(bbottom_node)), real(x(bbottom_node)),'k.','linewidth',3)


[up, zp, xp] = create2darray(real(z), real(x), real(u), delh, delh);
figure, clf, whitebg('white')
imagesc(zp, xp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('Field (real)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(real(z(btop_node)), real(x(btop_node)),'k.','linewidth',3)
plot(real(z(bbottom_node)), real(x(bbottom_node)),'k.','linewidth',3)

figure, clf, whitebg('white')
contour(zp, xp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('z (m)'); ylabel('x (m)');
title('Field (real)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(real(z(btop_node)), real(x(btop_node)),'k.','linewidth',3)
plot(real(z(bbottom_node)), real(x(bbottom_node)),'k.','linewidth',3)

disp('THE END !');

toc
