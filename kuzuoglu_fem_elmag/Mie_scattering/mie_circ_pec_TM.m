function [Escat] = mie_circ_pec_TM(freq, a, phii, phio, r)
% ***********************************************************************
% Scattered far-field of a PEC circular object (Mie series solution)
% TM_z polarization
%
% INPUT:
% freq: frequency in MHz
% a   : radius of the cylinder (m)
% phii: Angle of incidence (radian)
% phio: Angle of observation point (radian)
% r   : Radius of observation point (distance) (m)
% OUTPUT:
% Escat: scattered field
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

E0 = besselj(0,ka)*besselh(0,2,kr)/besselh(0,2,ka);
n = 1:50; % terms in the summation
E = sum(2*((-1j).^n).*cos(n*(phio-phii)).*besselj(n,ka).*besselh(n,2,kr)./besselh(n,2,ka));

Escat = -(E0+E); % scattered field
