function [qual] = triangle_quality1(x,y)
% It computes the quality measure Q1 for triangular elements.
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

sum1 = l12+l23+l31;
sum2 = -l12+l23+l31;
sum3 = l12-l23+l31;
sum4 = l12+l23-l31;

% radius of circumscribed circle
rcc = l12.*l23.*l31 ./ sqrt(sum1.*sum2.*sum3.*sum4);

% radius of inscribed circle
ric = 0.5*sqrt(sum2.*sum3.*sum4./sum1);

% quality measure
qual = 2*ric./rcc;
