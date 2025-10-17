% Name: FEM2D_surf.m 
% 2D Scattering from a Rough Surface (Sec. 5.6.6)
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
%***********
c0  = 3*1e8;           % m/sec, velocity of light in free space
nu0 = 120*pi;          % ohm, intrinsic impedance of the free space
e0  = (1e-9)/(36*pi);  % F/m, permittivity of free space
mu0 = 4*pi*1e-7;       % H/m, permeability of free space

% Input parameters
%******************
freq    = 300;          % MHz, frequency
freq    = freq*1e6;     % Hz, frequency
lambda0 = c0/freq;      % meter, wavelength
k0      = 2*pi/lambda0; % 1/meter, wavenumber
omg     = 2*pi*freq;    % rad/sec, radial frequency 

% Choose object type
surf_type = 'imp'; % 'pec' or 'imp' (PEC or impedance surface)

% Choose polarization
polrz = 'tm'; % 'TM' or 'TE' polarization

phii = 30;           % degree, angle of incident field
phii = phii*pi/180;  % radian, angle of incident field

% used for impedance surface
er = 80;     % dielectric constant
sigma = 4.8; % conductivity

% Mesh Generation
%******************
delh = lambda0/10;   % element size

surface_model = 'sea';  % 'flat', 'sea' or 'ground'

xdiscard = 1.5*lambda0; % discarded region from the left and right hand 
                        % side of the surface (used for postprocessing)
Len = 20*lambda0 + 2*xdiscard; % length of the rough surface

v = 10; % wind velocity (m/s) (used for sea surface)

h = .3*lambda0;    % rms height (used for ground surface)
lc = .3*lambda0;   % correlation length (used for ground surface)

obj = struct('delh',delh,'surface_model',surface_model,'Len',Len, ...
             'v',v,'h',h,'lc',lc, 'xdiscard', xdiscard); 

disp('The mesh is being created...');

[conn,x,y,pmlbin_node,pmlbout_node,pml_node,scatb_node,scatb_elm,huygb_node,huygb_elm,wind] ...
      = meshgen_surf(obj,freq);

disp('The mesh has been created!');

M = size(conn,1); % number of elements
N = length(x);    % number of nodes

% Material and source parameters
%*******************************
epsr = ones(M,1);  % relative permittivity of elements
mur = 1; % relative permeability

% used for impedance surface
epsr1 = er - 1j*60*sigma*lambda0;
if strcmpi(polrz, 'tm') % TM
    alfa2 = -1j*k0*(epsr1-1)^(1/2);
else
    alfa2 = -1j*k0*sqrt(epsr1-1)/epsr1;
end

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
disp('The matrix is being constructed...');
Ne = 3;           % number of nodes in each triangular element
A = sparse(N,N);  % initialize the global matrix to sparse zero matrix (size: NxN)
b = sparse(N,1);  % initialize the global b vector to sparse zero vector (size: Nx1)

% Start the matrix assembly process:
for e = 1:M  % loop over all elements
    if strcmpi(polrz, 'tm') % TM
        pe = 1/mur;
        qe = -k0^2*epsr(e);
    else % TE
        pe = 1/epsr(e);
        qe = -k0^2*mur;
    end

    %[Ae, ~, ~] = element_matrix_tri2(e, x, y, conn, pe, pe, qe, 0);    % element matrix and b-vector
    [Ae, ~, ~] = element_matrix_tri(e, x, y, conn, pe, pe, qe, 0);    % element matrix and b-vector

    for i = 1:Ne            % loop over the local nodes of each element
        ig = conn(e,i);     % global node corresponding to i
        for j = 1:Ne        % loop over the local nodes of each element
            jg = conn(e,j); % global node corresponding to j
            A(ig,jg) = A(ig,jg) + Ae(i,j); % fill the global matrix
        end
    end
end

disp('The matrix has been constructed!');

xr = real(x); 
yr = real(y); 

uinc = exp(1j*k0*(xr*cos(phii)+yr*sin(phii))); % incident field at nodes

% Imposition of BCs
%**********************
if strcmpi(surf_type, 'pec') % PEC object
    disp('The BCs are being imposed...');

    if strcmpi(polrz, 'tm') % TMz
        BCnodes = [scatb_node pmlbout_node];  % boundary node IDs
        BCvalues = [-uinc(scatb_node) zeros(size(pmlbout_node))];  % function values at the boundary nodes        
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
            
            [anx, any] = unit_vector(xr, yr, node1, node2, nodeout);
            anx = -anx; any = -any;
            
            delc = sqrt((xr(node1)-xr(node2))^2+(yr(node1)-yr(node2))^2);
            ym = 0.5*(yr(node1)+yr(node2));
            xm = 0.5*(xr(node1)+xr(node2));
            
            um = exp(1j*k0*(xm*cos(phii)+ym*sin(phii)));
            g = -1j*k0*um*(cos(phii)*anx + sin(phii)*any);
            term = g*delc*0.5;
            
            b(node1) = b(node1)+term;
            b(node2) = b(node2)+term;
            
        end
    end    
    disp('The BCs have been imposed!');
    
else
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
        
        [anx, any] = unit_vector(xr, yr, node1, node2, nodeout);
        anx = -anx; any = -any;
        
        delc = sqrt((xr(node1)-xr(node2))^2+(yr(node1)-yr(node2))^2);
        ym = 0.5*(yr(node1)+yr(node2));
        xm = 0.5*(xr(node1)+xr(node2));
        
        um = exp(1j*k0*(xm*cos(phii)+ym*sin(phii)));
        g = -1j*k0*um*(cos(phii)*anx + sin(phii)*any)  - alfa2*um;
        term = g*delc*0.5;
        
        b(node1) = b(node1)+term;
        b(node2) = b(node2)+term;
        
        A(node1,node1) = A(node1,node1) + (1/3)*alfa2*delc;
        A(node2,node2) = A(node2,node2) + (1/3)*alfa2*delc;
        A(node1,node2) = A(node1,node2) + (1/6)*alfa2*delc;
        A(node2,node1) = A(node2,node1) + (1/6)*alfa2*delc;
    end
end

b(scatb_node) = b(scatb_node).*wind';

 
% Solution of the global matrix equation
%***************************************
disp('The matrix equation is being solved...');

u = A\b;  % sattered field values at the nodes
u = u.';

utot = u + uinc; % total field
utot(pml_node) = 0;

disp('The matrix equation has been solved!');


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************

% RCS calculation
%****************
disp('The RCS is being calculated...');
 
dx = abs(max(xr)-min(xr));
dy = abs(max(yr)-min(yr));
Dmax = sqrt(dx^2+dy^2);     
rfar = 1000*(2*Dmax^2/lambda0); % far-field distance

if strcmpi(polrz, 'tm') % TM
    [RCS, farfield] = rcsTM(utot, conn, xr, yr, huygb_node, huygb_elm, 'diel');
else % TE
    [RCS, farfield] = rcsTE(utot, conn, xr, yr, huygb_node, huygb_elm, 'diel');
end

RCSdB = 10*log10(RCS/lambda0);

disp('The RCS has been calculated!');

% Plot the RCS
%*************
figure, clf, whitebg('white'), hold on
set(gcf, 'Color', [1 1 1])
plot(0:180, RCSdB(1:181),'k', 'Linewidth', [2]); axis tight;
title('Bistatic RCS Profile');
xlabel('\phi (degree)');
ylabel('RCS / \lambda (dB)');

% Plot the scattered field
%*************************
[up, xp, yp] = create2darray(xr, yr, real(u), delh, delh);

figure, clf, whitebg('white')
imagesc(xp, yp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Scattered Field (real)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(xr(scatb_node), yr(scatb_node),'k.','linewidth',3)

% Plot the total field
%*********************
[up, xp, yp] = create2darray(xr, yr, real(utot), delh, delh);

figure, clf, whitebg('white')
imagesc(xp, yp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Total Field (real)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(xr(scatb_node), yr(scatb_node),'k.','linewidth',3)

% Plot the incident field
%*********************
[up, xp, yp] = create2darray(xr, yr, real(uinc), delh, delh);

figure, clf, whitebg('white')
imagesc(xp, yp, up);
colormap(jet)
axis equal tight;
shading interp;
set(gca,'Ydir','normal');
colorbar;
xlabel('x (m)'); ylabel('y (m)');
title('Incident Field (real)');
set(gcf, 'Color', [1 1 1]);
hold on
plot(xr(scatb_node), yr(scatb_node),'k.','linewidth',3)

disp('THE END !');

toc
