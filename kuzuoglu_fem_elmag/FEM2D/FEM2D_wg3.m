%% Name: FEM2D_wg3.m 
% Line-source radiation in a parallel-plate waveguide (Sec. 5.6.9)
%
% The waveguide may or may not include an object.
% It solves the 2D Helmholtz equation with linear triangular elements.
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

% Constants
%**********
c0  = 3*1e8;           % m/sec, velocity of light in free space
nu0 = 120*pi;          % ohm, intrinsic impedance of the free space
e0  = (1e-9)/(36*pi);  % F/m, permittivity of free space
mu0 = 4*pi*1e-7;       % H/m, permeability of free space

% Input parameters
%*****************
freq    = 5000;         % MHz, frequency
freq    = freq*1e6;     % Hz, frequency
lambda0 = c0/freq;      % meter, wavelength
k0      = 2*pi/lambda0; % 1/meter, wavenumber
omg     = 2*pi*freq;    % rad/sec, radial frequency 

polrz = 'tm'; % 'TM' or 'TE' polarization

er = 4;  % dielectric constant of dielectric object (used if scat_type='diel')
mur = 1; % relative permittivity of object (used if scat_type='diel')

Is = 1e-4; % current on line source

% Mesh Generation
%****************
delh = 0.1/40;     % m, element size (lambda0/40)
dz = 0.5;          % m, length of the waveguide along z
dx = 0.1;          % m, length of the waveguide along x
object_exist = 0;  % flag showing an object exists or not within the waveguide
scat_type = 'pec'; % 'pec' or 'diel' (PEC or dielectric object)
scat_shape = 2;    % 1: circular, 2: rectangular
rc = 0.01;         % m, radius of circular object
lz = 0.1/2;        % m, half-length of rectangular object along z
lx = 0.1/4;        % m, half-length of rectangular object along x
cx = lx/2;         % m, center of the object along x
srcz = 0.05;       % m, line source location along z
srcx = dx/2;       % m, line source location along x
ls = 1*delh;       % m, half-length of source region
pmldim = 0.05;     % m, distance btw inner and outer PML boundaries

obj = struct('dz', dz, 'dx', dx, 'object_exist', object_exist, ...
    'scat_type', scat_type, 'scat_shape', scat_shape, 'rc', rc, ...
    'lz', lz, 'lx', lx, 'cx', cx, 'srcz', srcz, 'srcx', srcx, 'ls',ls, ...
    'pmldim', pmldim, 'delh', delh);

disp('The mesh is being created...');
[conn,z,x,pmlbinr_node,pmlboutr_node,pmlr_node, ...
    pmlbinl_node,pmlboutl_node,pmll_node,scatb_node,sourceb_node,...
    btop_node,bbottom_node,scatin_node,scatin_elm,sourceb_elm] ...
    = meshgen_wg3(obj);
disp('The mesh has been created!');
   
M = size(conn,1); % number of elements
N = length(z);    % number of nodes

% Cutoff frequency and wavenumbers
%*********************************
nmode = 1;    % mode
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

rho = sqrt((z(sourceb_node)-srcz).^2 + (x(sourceb_node)-srcx).^2);
term = -1.0*k0*nu0*0.25*Is*besselh(0, 2, k0*rho);  
     
if strcmpi(polrz, 'tm') % TMz
    if strcmpi(scat_type, 'pec')
       BCnodes = [scatb_node; sourceb_node; btop_node; bbottom_node];  % boundary node IDs
       BCvalues = [zeros(size(scatb_node)); term; zeros(size(btop_node)); zeros(size(bbottom_node))];  % function values at the boundary nodes        
    else
       BCnodes = [sourceb_node; btop_node; bbottom_node];  % boundary node IDs
       BCvalues = [term; zeros(size(btop_node)); zeros(size(bbottom_node))];  % function values at the boundary nodes
    end
else
    BCnodes = sourceb_node;  % boundary node IDs
    BCvalues = term;  % function values at the boundary nodes
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
title('Field (real)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(z(scatb_node), x(scatb_node),'k.','linewidth',3)

disp('THE END !');

toc
