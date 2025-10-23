function [Hscat] = mie_circ_pec_TE(freq, a, phii, phio, r)
% ***********************************************************************
% Scattered far-field of a PEC circular object (Mie series solution)
% TE_z polarization
%
% INPUT:
% freq: frequency in MHz
% a   : radius of the cylinder (m)
% phii: Angle of incidence (radian)
% phio: Angle of observation point (radian)
% r   : Radius of observation point (distance) (m)
% OUTPUT:
% Hscat: scattered field
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

freq   = freq*1e6;    % Hz, frequency
c0     = 3*1e8;       % m/sec, velocity of light in free space
lambda = c0/freq;     % meter, wavelength
k      = 2*pi/lambda; % 1/meter, wavenumber

phii = phii-pi;  

ka = k*a;
kr = k*r;

n = 0;
dbesselka = (besselj(n-1,ka)-besselj(n+1,ka))/2;  
dhankelka = (besselh(n-1,2,ka)-besselh(n+1,2,ka))/2;
H0 = dbesselka*besselh(0,2,kr)/dhankelka;

n = 1:50; % terms in the summation
dbesselka = (besselj(n-1,ka)-besselj(n+1,ka))/2;
dhankelka = (besselh(n-1,2,ka)-besselh(n+1,2,ka))/2;
H = sum(2*((-1j).^n).*cos(n*(phio-phii)).*dbesselka.*besselh(n,2,kr)./dhankelka);

Hscat = -(H0+H); % scattered field
