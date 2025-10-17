function [qual] = tetra_quality2(x,y,z)
% It computes the quality measure Q2 for tetrahedral elements.
% Input:
% x,y,z: vertex coordinates (each size: Nx4, where N is the # of tetrahedra)
% Output:
% qual: quality measure between 0 and 1 (size: Nx1)
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

v41 = [x(:,1)-x(:,4) y(:,1)-y(:,4) z(:,1)-z(:,4)];
v42 = [x(:,2)-x(:,4) y(:,2)-y(:,4) z(:,2)-z(:,4)];
v43 = [x(:,3)-x(:,4) y(:,3)-y(:,4) z(:,3)-z(:,4)];

v41norm2 = v41(:,1).^2 + v41(:,2).^2 + v41(:,3).^2;
v42norm2 = v42(:,1).^2 + v42(:,2).^2 + v42(:,3).^2;
v43norm2 = v43(:,1).^2 + v43(:,2).^2 + v43(:,3).^2;

v = v41norm2.*cross(v42,v43,2) + v42norm2.*cross(v43,v41,2) + v43norm2.*cross(v41,v42,2);
Z = sqrt(v(:,1).^2 + v(:,2).^2 + v(:,3).^2);

% volume of tetrahedral elements
vol = zeros(length(x),1);
for i = 1:length(x)
    vol(i) = det([1 1 1 1; x(i,1) x(i,2) x(i,3) x(i,4); ...
             y(i,1) y(i,2) y(i,3) y(i,4); z(i,1) z(i,2) z(i,3) z(i,4)])/6;
end

% radius of circumscribed sphere
rcs = Z ./ (12*vol);

% volume of circumscribed sphere
volcs = (4/3)*pi*rcs.^3;

% quality measure
qual = ((9*pi/(2*sqrt(3))) * (vol./volcs)).^(1/3);
