%% Name: FEM1D_wg.m 
% Propagation inside a Parallel-Plate Waveguide (Sec. 4.6.3)
%
% It solves 1-D Helmholtz’s equation with "linear" elements. 
% It considers both TM and TE modes. 
% The program displays the cutoff wavenumbers, cutoff frequencies and
% field solutions. Analytical values are also displayed and compared 
% with the numerical solution.
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
e0  = (1e-9)/(36*pi);  % F/m, permittivity of free space
mu0 = 4*pi*1e-7;       % H/m, permeability of free space

modeflag = 1; % if 1, TM mode, else TE mode

le = 1/40;    % m, element length
L = 1;        % m, total length of the domain (thickness of waveguide)

er = 2.7;     % dielectric constant of material within waveguide
mur = 1;      % relative permeability of material (non-magnetic medium)

% Mesh Generation
%****************
x = linspace(0, L, round(L/le)+1)';  % array of coordinates (size: Nx1)
xm = (x(1:end-1)+x(2:end))/2; % midpoint coords in the elements (size: Mx1)

N = length(x);   % number of nodes
M = N-1;         % number of linear elements
lev = diff(x);   % exact element length (size: Mx1)

conn = [(1:N-1)' (2:N)']; % connectivity matrix  (size: Mx2)
                           

%% ***********************************************************************
% MAIN BODY
%*************************************************************************
Ne = 2; % number of nodes in each element
A = sparse(N,N);  % initialize the global matrix to sparse zero matrix % (size: NxN)
B = sparse(N,N);  % initialize the global matrix to sparse zero matrix % (size: NxN)      

% Start the matrix assembly process:
%***********************************
for e = 1:M  % loop over all elements
    Ae = [1 -1;-1 1]/lev(e); % element matrix  
    Be = [2 1;1 2]*lev(e)/6; % element matrix      
    for i = 1:Ne            % loop over the local nodes of each element
        ig = conn(e,i);     % global node corresponding to i
        for j = 1:Ne        % loop over the local nodes of each element
            jg = conn(e,j); % global node corresponding to j
            A(ig,jg) = A(ig,jg) + Ae(i,j); % fill the global matrix
            B(ig,jg) = B(ig,jg) + Be(i,j); % fill the global matrix            
        end
    end
end

% Imposition of Dirichlet boundary condition in TM mode
%**********************************************************
if modeflag == 1  % TM
    A(1,:) = 0;   % make the first row of A matrix zero
    A(1,1) = 1;   % make the diagonal term of the A matrix 1
    A(N,:) = 0;   % make the last row of A matrix zero
    A(N,N) = 1;   % make the diagonal term of the A matrix 1
    
    B(1,:) = 0;   % make the first row of B matrix zero
    B(1,1) = 1;   % make the diagonal term of the B matrix 1
    B(N,:) = 0;   % make the last row of B matrix zero
    B(N,N) = 1;   % make the diagonal term of the B matrix 1    
    
    %A = A(1:end-1, 1:end-1);
    %A = A(2:end, 2:end);
    %B = B(1:end-1, 1:end-1);
    %B = B(2:end, 2:end);    
end

% Solution of the global matrix equation
%***************************************
[u,lambdae,iresult] = sptarn(A,B,-10,1500);

% eliminate meaningless values in the solutions
if modeflag == 1
   u = u(:, 3:end);
   lambdae = lambdae(3:end);
else
   u = u(:, length(find(lambdae < 1e-3))+1:end);
   lambdae = lambdae(lambdae > 1e-3);
end


%% ***********************************************************************
% POST-PROCESSING
%*************************************************************************
kc = real(sqrt(lambdae));             % cutoff wavenumber
fc = kc / (2*pi*sqrt(mur*mu0*er*e0)); % cutoff frequency
    
% Exact results
%***************
m = (1:length(kc))'; % vector of modes
kce = m*pi/L;        % exact cutoff wavenumber (rad/m)
fce = kce / (2*pi*sqrt(mur*mu0*er*e0)); % exact cutoff frequency (Hz)

% Percentage error
errk = 100*abs(kc-kce)./abs(kce); % error in kc
errf = 100*abs(fc-fce)./abs(fce); % error in fc

disp(' '); disp('Cutoff wavenumbers (rad/m)')
for i = 1:length(m)
   disp(['m = ' sprintf('%d\t', m(i)) 'kc = ' sprintf('%.3f\t', kc(m(i))) ...
         'kce = ' sprintf('%.3f\t', kce(m(i))) 'error (%) = ' sprintf('%.2f', errk(m(i)))]);
end

disp(' '); disp('Cutoff frequencies (Hz)')
for i = 1:length(m)
   disp(['m = ' sprintf('%d\t', m(i)) 'fc = ' sprintf('%.2e\t', fc(m(i))) ...
         'fce = ' sprintf('%.2e\t', fce(m(i))) 'error (%) = ' sprintf('%.2f', errf(m(i)))]);     
end


for i = 1:length(kce)
    if modeflag == 1 % TM
       ue(:,i) = sin(kce(i)*x); % exact field solution
    else %TE
       ue(:,i) = cos(kce(i)*x); % exact field solution        
    end
end

% Plot the results
%*****************
figure, clf, whitebg('white'), set(gcf, 'Color', [1 1 1]);
mp = 1:3;  % modes to be plotted
plot(x,u(:,mp) ./ max(max(u(:,mp))), 'Linewidth', 2) 
grid on
xlabel('x (m)')
if modeflag == 1
   ylabel('Normalized E_{0z} (V/m)')
else
   ylabel('Normalized H_{0z} (A/m)')    
end

figure, clf, whitebg('white'), set(gcf, 'Color', [1 1 1]);
mp = 1:3;  % modes to be plotted (exact field)
plot(x,ue(:,mp), 'Linewidth', 2) 
grid on
xlabel('x (m)')
if modeflag == 1
   ylabel('Normalized E_{0z} (V/m)')
else
   ylabel('Normalized H_{0z} (A/m)')    
end

% plot the error
figure, clf, whitebg('white'), set(gcf, 'Color', [1 1 1]);
plot(m, errk, 'Linewidth', 2) 
grid on
xlabel('m')
ylabel('Error in cutoff wavenumbers (%)')

toc
