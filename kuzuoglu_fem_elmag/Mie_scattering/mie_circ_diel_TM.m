function [Escat] = mie_circ_diel_TM(freq, a, epsr, phii, phio, r)
% ***********************************************************************
% Scattered far-field of a dielectric circular object (Mie series solution)
% TM_z polarization
%
% INPUT:
% freq: frequency in MHz
% a   : radius of the cylinder (m)
% epsr: dielectric constant of the cylinder
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

kd = sqrt(epsr)*k;
kda = kd*a;
ka = k*a;
kr = k*r;
    
nn = -50:50;
n = abs(nn);
dbesselkda = (besselj(n-1,kda)-besselj(n+1,kda))/2; 
dbesselka = (besselj(n-1,ka)-besselj(n+1,ka))/2;  
dhankelka = (besselh(n-1,2,ka)-besselh(n+1,2,ka))/2;
term1 = -besselj(n,ka)./besselh(n,2,ka);
term2 = epsr*dbesselkda./(kda*besselj(n,kda));
term3 = dbesselka./(ka*besselj(n,ka));
term4 = dhankelka./(ka*besselh(n,2,ka));
an = term1.*(term2-term3)./(term2-term4);

hankeltermn = sqrt(2*1j/(pi*kr))*exp(-1j*kr)*1j.^nn;
Escat = sum((1j.^(-nn)).*an.*hankeltermn.*exp(1j*nn*(phio-phii)));
