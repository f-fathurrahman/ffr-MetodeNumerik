function [Ae, be, vole] = element_matrix_tetra(e, x, y, z, conn, pxe, pye, pze, qe, fe)
% Element matrix equation with direct calculation for nodal tetrahedral
% elements 
% INPUT:
% e   : element ID
% x,y,z : x, y and z coordinates of all nodes in the mesh
% conn: nodal connectivity matrix
% pxe : px value of the e-th element
% pye : py value of the e-th element
% pze : pz value of the e-th element
% qe  : q value of the e-th element
% fe  : f value of the e-th element
% OUTPUT:
% Ae  : element matrix (4x4)
% be  : element right-hand-side vector (4x1)
% vole: volume of the e-th element
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

Ne = 4; % number of nodes in an element

n1 = conn(e,1); n2 = conn(e,2); n3 = conn(e,3); n4 = conn(e,4); % nodes of the element
x1 = x(n1); x2 = x(n2); x3 = x(n3); x4 = x(n4); % x-coordinates of the nodes
y1 = y(n1); y2 = y(n2); y3 = y(n3); y4 = y(n4); % y-coordinates of the nodes
z1 = z(n1); z2 = z(n2); z3 = z(n3); z4 = z(n4); % z-coordinates of the nodes

c1(1) = det([x2 x3 x4; y2 y3 y4; z2 z3 z4]);
c1(2) = -det([x1 x3 x4; y1 y3 y4; z1 z3 z4]);
c1(3) = det([x1 x2 x4; y1 y2 y4; z1 z2 z4]);
c1(4) = -det([x1 x2 x3; y1 y2 y3; z1 z2 z3]);

c2(1) = -det([1 1 1; y2 y3 y4; z2 z3 z4]);
c2(2) = det([1 1 1; y1 y3 y4; z1 z3 z4]);
c2(3) = -det([1 1 1; y1 y2 y4; z1 z2 z4]);
c2(4) = det([1 1 1; y1 y2 y3; z1 z2 z3]);

c3(1) = det([1 1 1; x2 x3 x4; z2 z3 z4]);
c3(2) = -det([1 1 1; x1 x3 x4; z1 z3 z4]);
c3(3) = det([1 1 1; x1 x2 x4; z1 z2 z4]);
c3(4) = -det([1 1 1; x1 x2 x3; z1 z2 z3]);

c4(1) = -det([1 1 1; x2 x3 x4; y2 y3 y4]);
c4(2) = det([1 1 1; x1 x3 x4; y1 y3 y4]);
c4(3) = -det([1 1 1; x1 x2 x4; y1 y2 y4]);
c4(4) = det([1 1 1; x1 x2 x3; y1 y2 y3]);

vole = sum(c1)/6; % volume of element

% Element matrix (first part)
% Calculate volume integral of delNi.delNj
Ae1 = zeros(Ne,Ne);
for i = 1:Ne
    for j = i:Ne
        Ae1(i,j) = pxe*c2(i)*c2(j) + pye*c3(i)*c3(j) + pze*c4(i)*c4(j);
        Ae1(j,i) = Ae1(i,j);
    end
end
Ae1 = Ae1 / (36*vole);

% Element matrix (second part)
Ae2 = qe*(vole/20)*[2 1 1 1; ...
                    1 2 1 1; ...
                    1 1 2 1; ...
                    1 1 1 2];

% Element matrix
Ae = Ae1 + Ae2;

% Element right-hand-side vector
be = [1;1;1;1]*fe*vole/4;
