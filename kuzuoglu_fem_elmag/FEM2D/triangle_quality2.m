function [qual] = triangle_quality2(x,y)
% It computes the quality measure Q2 for triangular elements.
% Input:
% x,y: vertex coordinates (each size: Nx3, where N is the # of triangles)
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

% edge lengths
l12 = sqrt((x(:,1)-x(:,2)).^2 + (y(:,1)-y(:,2)).^2);
l23 = sqrt((x(:,2)-x(:,3)).^2 + (y(:,2)-y(:,3)).^2);
l31 = sqrt((x(:,3)-x(:,1)).^2 + (y(:,3)-y(:,1)).^2);

% area
area = 0.5*abs((x(:,2)-x(:,1)).*(y(:,3)-y(:,1))-(x(:,3)-x(:,1)).*(y(:,2)-y(:,1)));

% quality measure
qual = 4*sqrt(3)*area./(l12.^2 + l23.^2 + l31.^2);