%% Name: FEM1D_slab.m 
% Scattering from a Lossy Dielectric Slab (Sec. 4.6.1)

% It solves 1-D Helmholtz’s equation with "linear" elements. 
% The slab is illuminated by a plane wave propagating along +x direction.
% Perpendicular polarization is considered.
% The code plots the scattered and total fields, and displays the
% magnitudes of reflection and transmission coefficients.
% Analytical values of the magnitudes of reflection and transmission 
% coefficients are also displayed and compared with the numerical solution.
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
c0  = 3*1e8;           % m/sec, velocity of light in free space
nu0 = 120*pi;          % ohm, intrinsic impedance of the free space
e0  = (1e-9)/(36*pi);  % F/m, permittivity of free space
mu0 = 4*pi*1e-7;       % H/m, permeability of free space

freq   = 300*1e6;      % Hz, frequency
lambda = c0/freq;      % meter, wavelength in free-space
k0     = 2*pi/lambda;  % 1/meter, wavenumber in free-space
omg    = 2*pi*freq;    % radial frequency 

le = lambda/40; % element length (in lambda) (initial)

L = 2;        % m, total length of the domain
d = 0.25;     % m, thickness of the slab

er = 2.7;     % dielectric constant (relative permittivity) of the slab
sigma = 5e-3; % S/m, conductivity of the slab
mur = 1;      % relative permeability of the slab (non-magnetic medium)

erc = er - 1j*sigma/(omg*e0); % complex dielectric constant of the slab

xa = L/2-d/2;  % x-coordinate of the slab (left)
xb = L/2+d/2;  % x-coordinate of the slab (right)

% Mesh Generation
%****************
x1 = linspace(0, xa, round(xa/le)+1)'; % coords in free-space on the left
x2 = linspace(xa, xb, round(d/le)+1)'; % coords in the slab
x3 = linspace(xb, L, round((L-xb)/le)+1)'; %coords in free-space on the right

x = [x1; x2(2:end); x3(2:end)]; % array of coordinates (size: Nx1)
xm = (x(1:end-1)+x(2:end))/2; % midpoint coords in the elements (size: Mx1)

N = length(x);   % number of nodes
M = N-1;         % number of linear elements
lev = diff(x);   % exact element length (size: Mx1)

conn = [(1:N-1)' (2:N)']; % connectivity matrix  (size: Mx2)
                                 
% Find the node and element IDs inside the slab
slabnodes = find((x > xa-1e-8) & (x < xb+1e-8)); % nodes inside the slab
slabelms = find((xm > xa) & (xm < xb));          % elements inside the slab

% Material and source parameters
%*******************************
evec = ones(M,1);      % initialization of vector of dielectric constant
evec(slabelms) = erc;  % assign erc values of the slab
muvec = ones(M,1);     % initialization of vector of relative permeability
muvec(slabelms) = mur; % assign mur values of the slab

uinc = exp(-1j*k0*xm); % incident electric field in each element (size:Mx1)

p = 1./muvec;          % p function in each element (size: Mx1)
q = -k0^2*evec;        % q function in each element (size: Mx1)
f = (-k0^2*p-q).*uinc; % f function in each element (size: Mx1)


%% ***********************************************************************
% MAIN BODY
%*************************************************************************
Ne = 2; % number of nodes in each element
A = sparse(N,N);  % initialize the global matrix to sparse zero matrix % (size: NxN)
b = sparse(N,1);  % initialize the global b vector to sparse zero vector % (size: Nx1)

% Start the matrix assembly process:
%***********************************
for e = 1:M  % loop over all elements
    Ae = [1 -1;-1 1]*p(e)/lev(e) + [2 1;1 2]*q(e)*lev(e)/6; % element matrix    
    be = [1;1]*0.5*lev(e)*f(e); % element b-vector
    for i = 1:Ne            % loop over the local nodes of each element
        ig = conn(e,i);     % global node corresponding to i
        for j = 1:Ne        % loop over the local nodes of each element
            jg = conn(e,j); % global node corresponding to j
            A(ig,jg) = A(ig,jg) + Ae(i,j); % fill the global matrix
        end
        b(ig) = b(ig) + be(i); % fill the global b vector
    end
end

% Imposition of absorbing boundary conditions
%*********************************************
A(1,1) = A(1,1) + 1j*k0*p(1);   % ABC at x=0
A(N,N) = A(N,N) + 1j*k0*p(end); % ABC at x=L

% Solution of the global matrix equation
%***************************************
u = A\b;     % scattered electric field values at the nodes (size: Nx1)
u = full(u); % convert the sparse array to full array


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************
uincn = exp(-1j*k0*x); % incident field at each node (size: Nx1)
utot = u + uincn;      % total field at each node (size: Nx1)   

% Magnitudes of reflection and transmission coefficients:
R = abs(u(1));      % reflection coefficient
T = abs(utot(end)); % transmission coefficient

disp(['|R| = ' sprintf('%.4f', R)]);
disp(['|T| = ' sprintf('%.4f', T)]);

% Exact reflection and transmission coefficients:
nu = nu0/sqrt(erc/mur); % ohm, intrinsic impedance of the slab
k = k0*sqrt(mur*erc);   % 1/meter, wavenumber in the slab

r = (nu-nu0)/(nu+nu0);
Re = r*(1-exp(-2j*k*d)) / (1-r^2*exp(-2j*k*d));
Te = (1-r^2)*exp(-1j*k*d) / (1-r^2*exp(-2j*k*d));
Re = abs(Re);  % exact reflection coefficient (mag)
Te = abs(Te);  % exact transmission coefficient (mag)

disp(['|Re| = ' sprintf('%.4f', Re)]);
disp(['|Te| = ' sprintf('%.4f', Te)]);

% Percentage error
errR = 100*abs(R-Re)/Re;
disp(['Error in R (%) = ' sprintf('%.2f', errR)]);

errT = 100*abs(T-Te)/Te;
disp(['Error in T (%) = ' sprintf('%.2f', errT)]);

% Plot the results
%*****************
figure, clf, whitebg('white'), set(gcf, 'Color', [1 1 1]);
fill([xa xb xb xa], [0 0 2 2], 'y'); hold on
plot(x,abs(u),'k', 'Linewidth', 2) % plot the scattered field
hold on
plot(x,abs(utot),'k-.', 'Linewidth', 2) % plot the total field
grid on
xlabel('x (m)')
ylabel('|Ez| (V/m)')
legend('Slab region', 'Scattered field', 'Total field')

toc