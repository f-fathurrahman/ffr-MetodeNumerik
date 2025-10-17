function [Ephi, Etheta] = mie_sphere_pec(freq, a, phio, thetao,  r)
% ***********************************************************************
% Scattered far-field of a PEC sphere (Mie series solution)
% The incident electric field is polarized along the x-direction and
% propagates along the +z-direction, i.e., Einc_x = exp(-jkz)
%
% INPUT:
% freq  : frequency in Hz
% a     : radius of the cylinder (m)
% phio  : Phi angle of observation point (radian)
% thetao: Theta angle of observation point (radian)
% r     : Radius of observation point (far-field distance) (m)
% OUTPUT:
% Ephi, Etheta: Phi and theta components of the scattered field
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

ka = k*a;

n = 1:50; % terms in the summation

dbesselka = (besselj(n-1/2,ka)-besselj(n+3/2,ka))/2;
besselka = besselj(n+1/2,ka);
dsbesselka = (pi/4)*((pi*ka/2)^(-1/2))*besselka + sqrt(pi*ka/2)*dbesselka;
sbesselka = sqrt(pi*ka/2)*besselka;

dhankelka = (besselh(n-1/2,2,ka)-besselh(n+3/2,2,ka))/2;
hankelka = besselh(n+1/2,2,ka);
dshankelka = (pi/4)*((pi*ka/2)^(-1/2))*hankelka + sqrt(pi*ka/2)*dhankelka;
shankelka = sqrt(pi*ka/2)*hankelka;

an = 1j.^(-n).*(2*n+1)./(n.*(n+1));
bn = -an.*dsbesselka./dshankelka;
cn = -an.*sbesselka./shankelka;

legn1 = zeros(1,length(n));
legn2 = zeros(1,length(n));
for ii = 1:length(n)
    legn = legendre(ii,cos(thetao));
    legn1(ii) = legn(2);
    if length(legn)<3
       legn2(ii) = 0;
    else
       legn2(ii) = legn(3);
    end
end
dlegn1 = (-1/sin(thetao))*(legn2+(cos(thetao)/sin(thetao))*legn1);

if thetao == pi
   term1 = ((-1).^n/2).*n.*(n+1);    
   term2 = 1j.^n.*(bn-cn).*term1;
   Etheta = sum(term2);
   Ephi = Etheta;               
else
   Etheta = sum(1j.^n.*(bn*sin(thetao).*dlegn1-cn.*legn1/sin(thetao)));
   Ephi = sum(1j.^n.*(bn.*legn1/sin(thetao)-cn.*dlegn1*sin(thetao)));
end

Etheta = Etheta*cos(phio)*1j*exp(-1j*k*r)/(k*r);
Ephi = Ephi*sin(phio)*1j*exp(-1j*k*r)/(k*r);
