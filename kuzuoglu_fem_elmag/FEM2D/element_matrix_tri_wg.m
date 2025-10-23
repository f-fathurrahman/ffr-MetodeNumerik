function [Ae1, Ae2] = element_matrix_tri_wg(e, x, y, conn, pxe, pye, qe)
% Element matrix equation for waveguide application (direct calculation) 
% (triangular elements)
% INPUT:
% e   : element ID
% x,y : x and y coordinates of all nodes in the mesh
% conn: connectivity matrix
% pxe : pxx value of the e-th element
% pye : pyy value of the e-th element
% qe  : q value of the e-th element
% OUTPUT:
% Ae1  : element matrix 1 (3x3)
% Ae2  : element matrix 2 (3x3)
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

n1 = conn(e,1);  n2 = conn(e,2);  n3 = conn(e,3);
x1 = x(n1);  x2 = x(n2);  x3 = x(n3);
y1 = y(n1);  y2 = y(n2);  y3 = y(n3);

delxdelpsi = x2-x1;
delxdelnu = x3-x1;
delydelpsi = y2-y1;
delydelnu = y3-y1;

Jdet = delxdelpsi*delydelnu-delxdelnu*delydelpsi; %determinant of jacobian
% areae = 0.5*Jdet; % area of the element

Ae1(1,1) = pxe*(y2-y3)^2 + pye*(x3-x2)^2;
Ae1(1,2) = pxe*(y2-y3)*(y3-y1) + pye*(x3-x2)*(x1-x3);
Ae1(1,3) = pxe*(y2-y3)*(y1-y2) + pye*(x3-x2)*(x2-x1);
Ae1(2,1) = Ae1(1,2);
Ae1(2,2) = pxe*(y3-y1)^2 + pye*(x1-x3)^2;
Ae1(2,3) = pxe*(y3-y1)*(y1-y2) + pye*(x1-x3)*(x2-x1);
Ae1(3,1) = Ae1(1,3);
Ae1(3,2) = Ae1(2,3);
Ae1(3,3) = pxe*(y1-y2)^2 + pye*(x2-x1)^2;
Ae1 = Ae1/(2*Jdet);

Ae2(1,1) = qe*Jdet/12;
Ae2(1,2) = qe*Jdet/24;
Ae2(1,3) = qe*Jdet/24;
Ae2(2,1) = Ae2(1,2);
Ae2(2,2) = qe*Jdet/12;
Ae2(2,3) = qe*Jdet/24;
Ae2(3,1) = Ae2(1,3);
Ae2(3,2) = Ae2(2,3);
Ae2(3,3) = qe*Jdet/12;

