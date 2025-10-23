function [Ae, be, areae] = element_matrix_quad(e, x, y, conn, pxe, pye, qe, fe)
% Element matrix equation with MATLAB's numerical integration function
% (quadrilateral elements) 
% INPUT:
% e    : element ID
% x,y  : x and y coordinates of all nodes in the mesh
% conn : connectivity matrix
% pxe  : px value of the e-th element
% pye  : py value of the e-th element
% qe   : q value of the e-th element
% fe   : f value of the e-th element
% OUTPUT:
% Ae   : element matrix (3x3)
% be   : element right-hand-side vector (3x1)
% areae: area of the e-th element
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

global connq;
global xq; global yq;
global i; global j;
global pxeq; global pyeq;
global qeq; global feq;
global ee;

connq = conn; ee = e;
xq = x; yq = y;
pxeq = pxe; pyeq = pye;
qeq = qe; feq = fe;

Ne = 4; % number of nodes in an element

Ae = zeros(Ne, Ne);         
be = zeros(Ne, 1);         

% Element matrix
for i = 1:Ne
    for j = i:Ne        
        Ae(i,j) = integral2(@integrnd_quad_Ae, -1, 1, -1, 1);
        Ae(j,i) = Ae(i,j);
    end
end

% Element right-hand-side vector
for i = 1:Ne
    be(i) = integral2(@integrnd_quad_be, -1, 1, -1, 1);
end

% Area of the element
areae = integral2(@integrnd_quad_areae, -1, 1, -1, 1);
