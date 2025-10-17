function [Ae,be,vole] = element_matrix_hexavec(e,x,y,z,Len,conn,...
                        elm2edge_conn, pxe,pye,pze,qxe,qye,qze,Fxe,Fye,Fze)
% Element matrix equation with Gaussian quadrature method for
% hexahedral edge (vector) elements 
% INPUT:
% e : element ID
% x,y,z : x, y and z coordinates of all nodes in the mesh
% Len: length of all edges in the mesh
% conn: element-to-node connectivity matrix
% elm2edge_conn: element-to-edge connectivity matrix
% edge2node_conn: edge-to-node connectivity matrix
% pxe : px value of the corresponding element
% pye : py value of the corresponding element
% pze : pz value of the corresponding element
% qxe : qx value of the corresponding element
% qye : qy value of the corresponding element
% qze : qz value of the corresponding element
% Fxe : Fx value of the corresponding element
% Fye : Fy value of the corresponding element
% Fze : Fz value of the corresponding element
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

global connq;
global xq; global yq; global zq;
global i; global j; global ee;
global pxeq; global pyeq; global pzeq;
global qxeq; global qyeq; global qzeq;
global Fxeq; global Fyeq; global Fzeq;
global elm2edge_connq; global Lenq;

connq = conn; ee = e;
xq = x; yq = y; zq = z;
pxeq = pxe; pyeq = pye; pzeq = pze;
qxeq = qxe; qyeq = qye; qzeq = qze;
Fxeq = Fxe; Fyeq = Fye; Fzeq = Fze;
elm2edge_connq = elm2edge_conn;
Lenq = Len;

Ne = 12; % number of edges in an element

Ae = zeros(Ne, Ne);         
be = zeros(Ne, 1);         

% Gauss weights and points (8-point)
gw = [1 1 1 1 1 1 1 1];   % Gauss weights 
a = 1/sqrt(3);
ksi = [-a -a -a -a a a a a];    % Gauss point ksi
nu = [-a -a a a -a -a a a];     % Gauss point nu
zeta = [-a a -a a -a a -a a];

% Element matrix
for i = 1:Ne
    for j = i:Ne
        sumg = 0;
        for k = 1:length(gw)                 
            sumg = sumg+gw(k)*integrnd_hexavec_Ae(ksi(k), nu(k), zeta(k));
        end       
        Ae(i,j) = sumg;
        Ae(j,i) = Ae(i,j);
    end
end

% Element right-hand-side vector
for i = 1:Ne
    sumg = 0;
    for k = 1:length(gw)       
        sumg = sumg+gw(k)*integrnd_hexavec_be(ksi(k), nu(k), zeta(k));
    end
    be(i) = sumg;
end

% Volume of the element
sumg = 0;
for k = 1:length(gw)
    sumg = sumg+gw(k)*integrnd_hexa_vole(ksi(k), nu(k), zeta(k));
end
vole = sumg;
