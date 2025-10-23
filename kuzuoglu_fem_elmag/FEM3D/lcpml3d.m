function [xc, yc, zc] = lcpml3d(x, y, z, k, pmlbin_x, pmlbin_y, pmlbin_z, ...
                                         pmlbout_x, pmlbout_y, pmlbout_z)
% It implements the Locally-Conformal PML method in 3D.
% INPUT:
% x,y,z: real coordinates (a single point)
% k: wavenumber
% pmlbin_x, pmlbin_y, pmlbin_z: coordinates of the points on the inner PML boundary
% pmlbout_x, pmlbout_y, pmlbout_z: coordinates of the points on the outer PML boundary
% OUTPUT:
% xc, yc, zc: complex coordinates (a single point)
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

% LC-PML parameters
alpha = 7.0*k;
alphajk = alpha/(1j*k);
m = 3;   % PML decay rate, integer 2 or 3
 
% find the point on the inner PML boundary which is nearest to point P
[ksi,ind] = min(sqrt((pmlbin_x-x).^2+(pmlbin_y-y).^2+(pmlbin_z-z).^2));
x0 = pmlbin_x(ind);
y0 = pmlbin_y(ind);
z0 = pmlbin_z(ind);

% find the point on the outer PML boundary in the direction of the unit vector
vpx = pmlbout_x-x0;
vpy = pmlbout_y-y0;
vpz = pmlbout_z-z0;
lp = sqrt(vpx.^2 + vpy.^2 + vpz.^2);
npx = vpx ./ lp; % unit vector from r0 to r1 (x-comp)
npy = vpy ./ lp; % unit vector from r0 to r1 (y-comp)
npz = vpz ./ lp; % unit vector from r0 to r1 (z-comp)

vx = x-x0;
vy = y-y0;
vz = z-z0;
l = sqrt(vx.^2 + vy.^2 + vz.^2);
nx = vx / l; % unit vector from r0 to r (x-comp)
ny = vy / l; % unit vector from r0 to r (y-comp)
nz = vz / l; % unit vector from r0 to r (z-comp)

if l<1e-8
    xc = x; yc = y; zc = z;
else
    [~,ind] = min(acos(nx*npx+ny*npy+nz*npz));
    x1 = pmlbout_x(ind);
    y1 = pmlbout_y(ind);
    z1 = pmlbout_z(ind);
    
    % local PML thickness, distance btw P0 and P1
    dpml = sqrt((x1-x0)^2+(y1-y0)^2+(z1-z0)^2);
    
    term = alphajk*((ksi^m) / (m*(dpml^(m-1))));
    
    % complex coordinates
    xc = x + term*nx;
    yc = y + term*ny;
    zc = z + term*nz;    
end
