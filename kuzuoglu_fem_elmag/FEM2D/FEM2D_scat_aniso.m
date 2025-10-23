%% Name: FEM2D_scat_aniso.m 
% 2D Scattering from an Anisotropic Object (Sec. 5.6.4)
%
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

global e0;
global nu0;
global mu0;
global k0;
global omg;
global rfar;

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

% Choose polarization
polrz = 'tm'; % 'TM' or 'TE' polarization

mur = 1;   % relative permittivity of anisotropic object
erxx = 9;  % xx comp. of dielectric constant of anisotropic object
eryy = 4;  % yy comp. of dielectric constant of anisotropic object
erzz = 2;  % zzz comp. of dielectric constant of anisotropic object
erxy = 0;  % xy comp. of dielectric constant of anisotropic object
eryx = 0; % yx comp. of dielectric constant of anisotropic object

er = [erxx erxy 0; eryx eryy 0; 0 0 erzz];

phii = 0;            % degree, angle of incident field
phii = phii*pi/180;  % radian, angle of incident field

scat_type = 'diel'; % do NOT change the type !!!

% Mesh Generation
%****************
delh = lambda0/40;   % element size

scat_shape = 1; 
% 1: circular object with circular domain
% 2: circular object with rectangular domain
% 3: half-circular object with rectangular domain
% 4: elliptical object with rectangular domain 
% 5: polygonal object with rectangular domain 
% 6: conecircle object with rectangular domain 
% 7: ogive object with rectangular domain 

disp('The mesh is being created...');
if (scat_shape <= 3)
    rc = 1*lambda0;      % radius of circular or half-circular scatterer
    obj = struct('rc', rc,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape);  
 
elseif scat_shape == 4
    rex = 1*lambda0;      % x-dim of ellipse
    rey = 0.5*lambda0;    % y-dim of ellipse    
    obj = struct('rex', rex, 'rey', rey,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape);  
   
elseif scat_shape == 5
    svx = [1.5 1 1 -1 -1 1.5];   % Scatterer Vertex Coordinates (x)
    svy = [1 1 -1 -1 -1.5 -1.5]; % Scatterer Vertex Coordinates (y)
    obj = struct('svx', svx, 'svy', svy,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape);    
   % Vertices must be in the counter-clockwise direction!
        
elseif scat_shape == 6
    rc = 1.6*lambda0;     % radius of circular base
    hang = 30 *pi/180;    % half cone angle
    hh = rc/tan(hang);    % height of the cone
    obj = struct('rc', rc, 'hang', hang, 'hh', hh,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape);    

elseif scat_shape == 7
    hang = 25 *pi/180;          % angle
    hh = 2.0*lambda0;             % half height
    rc = hh / (1/tan(hang*0.5));  % half height
    obj = struct('rc', rc, 'hang', hang, 'hh', hh,'delh',delh,'scat_type', scat_type, 'scat_shape', scat_shape);    
end

if scat_shape == 1
    [conn,x,y,pmlbin_node,pmlbout_node,pml_node,scatb_node, ...
        scatin_node,scatin_elm,scatb_elm,huygb_node,huygb_elm] ...
        = meshgen_scat_circle(obj,freq);    
else
    [conn,x,y,pmlbin_node,pmlbout_node,pml_node,...
        scatb_node,scatin_node,scatin_elm,scatb_elm,huygb_node,huygb_elm] ...
        = meshgen_scat(obj,freq);
end
  
disp('The mesh has been created!');

M = size(conn,1); % number of elements
N = length(x);    % number of nodes

% Material and source parameters
%*******************************
epsrxx = ones(M,1);    % relative permittivity of elements
epsryy = ones(M,1);    % relative permittivity of elements
epsrzz = ones(M,1);    % relative permittivity of elements
epsrxy = zeros(M,1);   % relative permittivity of elements
epsryx = zeros(M,1);   % relative permittivity of elements
epsrxx(scatin_elm) = erxx;  
epsryy(scatin_elm) = eryy;  
epsrzz(scatin_elm) = erzz;  
epsrxy(scatin_elm) = erxy;  
epsryx(scatin_elm) = eryx;  


%% ***********************************************************************
% MAIN BODY
%*************************************************************************

% Implement Locally-Conformal PML (LC-PML)
%*****************************************
pmlbin_x = x(pmlbin_node);    pmlbin_y = y(pmlbin_node);
pmlbout_x = x(pmlbout_node);  pmlbout_y = y(pmlbout_node);

disp('The LC-PML is being implemented...');
for i = 1:length(pml_node)
    in = pml_node(i);
    [x(in), y(in)] = lcpml(x(in), y(in), k0, pmlbin_x, pmlbin_y, pmlbout_x, pmlbout_y);
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
        pxxe = mur;
        pyye = pxxe;
        pxye = 0; pyxe = 0;
        qe = -k0^2*epsrzz(e);
    else % TE
        pxxe = epsrxx(e);
        pyye = epsryy(e);
        pxye = epsrxy(e);
        pyxe = epsryx(e);
        qe = -k0^2*mur;
    end

    [Ae, ~, ~] = element_matrix_tri2_aniso(e, x, y, conn, pxxe, pyye, pxye, pyxe, qe, 0); % element matrix and b-vector

    for i = 1:Ne            % loop over the local nodes of each element
        ig = conn(e,i);     % global node corresponding to i
        for j = 1:Ne        % loop over the local nodes of each element
            jg = conn(e,j); % global node corresponding to j
            A(ig,jg) = A(ig,jg) + Ae(i,j); % fill the global matrix
        end
    end
end

uinc = exp(1j*k0*(real(x)*cos(phii)+real(y)*sin(phii))); % incident field at nodes

Au = A*uinc;
b(scatin_node) = -1.0*Au(scatin_node);

disp('The matrix has been constructed!');

% Solution of the global matrix equation
%***************************************
disp('The matrix equation is being solved...');

u = A\b; % sattered field values at the nodes

utot = u + uinc; % total field
utot(pml_node) = 0;

disp('The matrix equation has been solved!');


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************

% RCS calculation
%****************
disp('The RCS is being calculated...');
 
dx = abs(max(x)-min(x));
dy = abs(max(y)-min(y));
Dmax = sqrt(dx^2+dy^2);     
rfar = 1000*(2*Dmax^2/lambda0); % far-field distance

if strcmpi(polrz, 'tm') % TM
  [RCS, farfield] = rcsTM(utot, conn, x, y, huygb_node, huygb_elm, scat_type);
else % TE
  [RCS, farfield] = rcsTE(utot, conn, x, y, huygb_node, huygb_elm, scat_type);   
end

RCSdB = 10*log10(RCS/lambda0);

disp('The RCS has been calculated!');

% Plot the RCS
%*************
figure, clf, whitebg('white'), hold on
set(gcf, 'Color', [1 1 1])
plot(0:359, RCSdB,'k', 'Linewidth', [2]); axis tight;
title('Bistatic RCS Profile');
xlabel('\phi (degree)');
ylabel('RCS / \lambda (dB)');

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
plot(x(scatb_node), y(scatb_node),'k.','linewidth',3)

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
plot(x(scatb_node), y(scatb_node),'k.','linewidth',3)

disp('THE END !');

toc
