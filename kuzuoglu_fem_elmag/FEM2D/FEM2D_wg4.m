%% Name: FEM2D_wg4.m 
% Propagation inside a 3-D waveguide with uniform cross section (Sec. 5.6.10)
%
% The cross-section of the waveguide can be rectangular or circular.
% It solves the 2D Helmholtz equation with linear triangular elements.
% It considers both TM and TE modes. 
% The program displays the cutoff wavenumbers, cutoff frequencies and the
% field solutions.
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
e0  = (1e-9)/(36*pi);  % F/m, permittivity of free space
mu0 = 4*pi*1e-7;       % H/m, permeability of free space

polrz = 'tm'; % 'TM' or 'TE'polarization

% Mesh Generation
%****************
delh = .5/40;  % m, element size
wg_shape = 2; % 1: circular, 2: rectangular
rc = 1;       % m, radius of circular waveguide
lx = 1;     % m, half-length of rectangular wavegide along x
ly = 0.5;     % m, half-length of rectangular wavegide along y

obj = struct('wg_shape', wg_shape, 'rc', rc, 'lx',lx,'ly',ly,'delh',delh);

disp('The mesh is being created...');
[conn,x,y,wgb_node] = meshgen_wg4(obj);
disp('The mesh has been created!');

M = size(conn,1); % number of elements
N = length(x);    % number of nodes


%% ***********************************************************************
% MAIN BODY
%*************************************************************************
Ne = 3;           % number of nodes in each triangular element
A = sparse(N,N);  % initialize the global matrix to sparse zero matrix (size: NxN)
B = sparse(N,N);  % initialize the global matrix to sparse zero matrix (size: NxN)

disp('The matrix is being constructed...');
% Start the matrix assembly process:
for e = 1:M  % loop over all elements
    %[Ae, Be] = element_matrix_tri2_wg(e, x, y, conn, 1,1,1);  % element matrix and b-vector
    [Ae, Be] = element_matrix_tri_wg(e, x, y, conn, 1,1,1);    % element matrix and b-vector
    
    for i = 1:Ne            % loop over the local nodes of each element
        ig = conn(e,i);     % global node corresponding to i
        for j = 1:Ne        % loop over the local nodes of each element
            jg = conn(e,j); % global node corresponding to j
            A(ig,jg) = A(ig,jg) + Ae(i,j); % fill the global matrix
            B(ig,jg) = B(ig,jg) + Be(i,j); % fill the global matrix            
        end
    end
end

disp('The matrix has been constructed!');

% Imposition of BCs 
%**********************
if strcmpi(polrz, 'tm') % TMz
    disp('The BCs are being imposed...');
    BCnodes = wgb_node;  % boundary node IDs
    BCvalues = zeros(size(wgb_node));  % function values at the boundary nodes
    
    % the following is a fast implementation of Dirichlet BCs
    nodes = 1:N;
    nodes(BCnodes) = 0;
    innernodes = find(nodes ~= 0);
    Atemp = speye(N, N);
    Atemp(innernodes,:) = 0;
    A(BCnodes,:) = 0;
    A = A+Atemp;
    B(BCnodes,:) = 0;
    B = B+Atemp;
    
    A = A(innernodes,:);
    A = A(:,innernodes);
    B = B(innernodes,:);
    B = B(:,innernodes);
    
    disp('The BCs have been imposed!');
else
    innernodes = 1:N;
end

% Solution of the global matrix equation
%***************************************
disp('The matrix equation is being solved...');

[u,lambdae,~] = sptarn(A,B,-10,200);

kc = real(sqrt(lambdae));      % cutoff wavenumber
fc = kc / (2*pi*sqrt(mu0*e0)); % cutoff frequency

if strcmpi(polrz, 'te') 
    kc = kc(2:end);
    fc = kc(2:end);
    u = u(:,2:end);  
end

disp('Cutoff wavenumbers:')
for m = 1:5
    for n = 1:5        
       disp(['kc = ' sprintf('%.4f', kc(m*n))]);
    end
end

disp('The matrix equation has been solved!');

% Plot the field
%***************
mp = 1:7;  % modes to be plotted

for i = 1:length(mp)
    [up, xp, yp] = create2darray(real(x(innernodes)), real(y(innernodes)), ...
                                 abs(u(:,mp(i))), delh, delh);
    figure, clf, whitebg('white'), set(gcf, 'Color', [1 1 1]);
    imagesc(xp, yp, up ./ max(max(up)));
    %contour(xp, yp, up ./ max(max(up)));    
    colormap(jet)
    axis equal tight;
    shading interp;
    set(gca,'Ydir','normal');
    colorbar;
    xlabel('x (m)'); ylabel('y (m)');
    if strcmpi(polrz, 'tm')
        title('Normalized E_{0z}')
    else
        title('Normalized H_{0z}')
    end
    hold on
    plot(x(wgb_node),y(wgb_node),'k.')
end

disp('THE END !');

toc
