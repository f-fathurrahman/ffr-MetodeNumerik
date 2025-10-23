function [xc, yc] = lcpml(x, y, k, pmlbin_x, pmlbin_y, pmlbout_x, pmlbout_y)
% It implements the Locally-Conformal PML method.
% INPUT:
% x,y: real coordinates (a single point)
% k: wavenumber
% pmlbin_x, pmlbin_y: coordinates of the points on the inner PML boundary
% pmlbout_x, pmlbout_y: coordinates of the points on the outer PML boundary
% OUTPUT:
% xc, yc: complex coordinates (a single point)
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
[ksi,ind] = min(sqrt((pmlbin_x-x).^2+(pmlbin_y-y).^2));
x0 = pmlbin_x(ind);
y0 = pmlbin_y(ind);
 
% find the point on the outer PML boundary in the direction of the unit vector
vpx = pmlbout_x-x0;
vpy = pmlbout_y-y0;
lp = sqrt(vpx.^2 + vpy.^2);
npx = vpx ./ lp; % unit vector from r0 to r1 (x-comp)
npy = vpy ./ lp; % unit vector from r0 to r1 (y-comp)

vx = x-x0;
vy = y-y0;
l = sqrt(vx.^2 + vy.^2);
nx = vx / l; % unit vector from r0 to r (x-comp)
ny = vy / l; % unit vector from r0 to r (y-comp)

if l<1e-8
    xc = x; yc = y;
else
    [~,ind] = min(acos(nx*npx+ny*npy));
    x1 = pmlbout_x(ind);
    y1 = pmlbout_y(ind);
    
    % local PML thickness, distance btw P0 and P1
    dpml = sqrt((x1-x0)^2+(y1-y0)^2);
    
    term = alphajk*((ksi^m) / (m*(dpml^(m-1))));
    
    % complex coordinates
    xc = x + term*nx;
    yc = y + term*ny;
end
