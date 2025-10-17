%% Name: FEM1D_app2.m 
% Application 2: Capacitor with Inhomogeneous Medium having 
% Spatially-Varying Permittivity Profile (Sec. 4.5.2)
%
% It solves 1-D generalized Laplace’s equation with "linear" elements 
% and Dirichlet boundary conditions.
% Permittivity profile is e0*(1+x/L).
% Analytical results are also computed.
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

% Input parameters and constants
%*******************************
L = 4e-2;      % total length of the domain (m)
M = 10;        % number of linear elements

e0 = 8.85e-12; % free-space permittivity (F/m)
rho0 = 50e-9;  % constant term in charge density (C/m)

V0 = 1;        % known potential at the first node (V)
VL = 0;        % known potential at the first node (V)

% Mesh Generation
%****************
N = M+1;      % number of nodes
le = L/M;     % element length
x = linspace(0, L, N)';          % array of coordinates (size: Nx1)
xm = linspace(le/2, L-le/2, M)'; % coordinates of the midpoints of the elements (size: Mx1)
conn = [(1:N-1)' (2:N)'];        % connectivity matrix  (size: Mx2)

% Material and source parameters
%*******************************
eps = e0*(1+xm/L);      % permittivity in each element (size: Mx1)
rhov = zeros(size(xm)); % zero array of charge density (size: Mx1)


%% ***********************************************************************
% MAIN BODY
%*************************************************************************
Ne = 2; % number of nodes in each element
A = sparse(N,N);  % initialize the global matrix to sparse zero matrix (size: NxN)
b = sparse(N,1);  % initialize the global b vector to sparse zero vector (size: Nx1)

% Start the assembly process:
%****************************
for e = 1:M  % loop over all elements
    be = [0;0]; % zero element b-vector
    %be = [1;1]*0.5*le*rhov(e);  % element b-vector (enable if rhov is nonzero)
    Ae = [1 -1; -1 1]*eps(e)/le; % element matrix    
    for i = 1:Ne            % loop over the local nodes of each element
        ig = conn(e,i);     % global node corresponding to i
        for j = 1:Ne        % loop over the local nodes of each element
            jg = conn(e,j); % global node corresponding to j
            A(ig,jg) = A(ig,jg) + Ae(i,j); % fill the global matrix
        end
        b(ig) = b(ig) + be(i); % fill the global b vector
    end
end

% Imposition of the Dirichlet boundary conditions (first approach)
%****************************************************************
BCnodes = [1 N];          % boundary node IDs
BCvalues = [V0 VL];       % function values at the boundary nodes
for i = 1:length(BCnodes) % loop over all boundary nodes
    n = BCnodes(i);       % boundary node 
    A(n,:) = 0;           % make the n-th row of A matrix zero
    A(n,n) = 1;           % make the diagonal term of the A matrix 1
    b(n) = BCvalues(i);   % assign the n-th row of b to the known value   
end

% Solution of the matrix system
%******************************
u = A\b;     % potential values at the nodes
u = full(u); % convert the sparse array to full array


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************

% Potential distribution
%***********************
num = 11; % number of ksi points within each element
ksi = linspace(-1,1,num)'; % array of chosen points in master domain
V = zeros(M*(num-1)+1, 1); % initialize the potential in entire domain
i = 1; j = 1;
for e = 1:M % loop over all elements
    V(i:i+num-1) = u(j)*0.5*(1-ksi) + u(j+1)*0.5*(1+ksi);
    i = i+num-1; 
    j = j+1;    
end
x2 = linspace(0,L,length(V))';   % coordinates in the new domain
xm2 = (x2(1:end-1)+x2(2:end))/2; % midpoints of x2

% Electric field within each element
%***********************************
Ex = -diff(u) / le;
 
% Exact potential and electric field in x2 and xm2
%*************************************************
Ve = V0+(VL-V0)*log(1+x2/L)/log(2);
Exe = (V0-VL)*(1./(xm2+L))/log(2);

% Error comparing the potentials found by FEM and the analytical method
%**********************************************************************
err = 100*norm(V-Ve).^2/norm(Ve).^2;
disp(['Error (%) = ' sprintf('%.4e', err)]);

% Error comparing the NODAL potentials found by FEM and the analytical method
%************************************************************************
ue = V0+(VL-V0)*log(1+x/L)/log(2);

err = 100*norm(u-ue).^2/norm(ue).^2;
disp(['Error (nodal values) (%) = ' sprintf('%.4e', err)]);

% Plot the results
%*****************
figure, clf, whitebg('white'), set(gcf, 'Color', [1 1 1]);
plot(x2,V,'k', 'Linewidth', 2)
hold on
plot(x2,Ve,'k--', 'Linewidth', 2)
grid on
xlabel('x (m)')
ylabel('Potential (V)')
legend('FEM', 'Analytical')
 
figure, clf, whitebg('white'), set(gcf, 'Color', [1 1 1]);
stairs(x,[Ex; Ex(end)],'k', 'Linewidth', 2)
hold on
plot(xm2,Exe,'k--', 'Linewidth', 2)
grid on
xlabel('x (m)')
ylabel('Electric field (V/m)')
legend('FEM', 'Analytical')

toc
