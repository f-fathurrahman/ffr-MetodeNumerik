function [Ae, be, vole] = element_matrix_tetravec(e, x, y, z, Len, conn, ...
                          elm2edge_conn, edge2node_conn, ...
                          pxe, pye, pze, qxe, qye, qze, Fxe, Fye, Fze)
% Element matrix equation with direct calculation for tetrahedral edge
% (vector) element 
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
% Ae  : element matrix (6x6)
% be  : element right-hand-side vector (6x1)
% vole: volume of the corresponding element
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

Ae1 = zeros(6,6);
Ae2 = zeros(6,6);
be = zeros(6,1);

elm_nodes = conn(e,:);
elm_edges = elm2edge_conn(e,:);
x = x(elm_nodes);
y = y(elm_nodes);
z = z(elm_nodes);
Len = Len(elm_edges);
edge2node_conne = edge2node_conn(elm_edges,:);

x1 = x(1); x2 = x(2); x3 = x(3); x4 = x(4); % x-coordinates of the nodes
y1 = y(1); y2 = y(2); y3 = y(3); y4 = y(4); % y-coordinates of the nodes
z1 = z(1); z2 = z(2); z3 = z(3); z4 = z(4); % z-coordinates of the nodes

local_edge_def = fun_local_edge_def(elm_nodes, edge2node_conne);

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

% function declarations
NiNj = @(i,j) (i==j).*(1/10) +(i~=j).*(1/20);

% Calculate Ae1: volume integral of pe.(delxNi).(delxNj)
% Calculate Ae2: volume integral of qe.Ni.Nj
for i = 1:6
    i1 = local_edge_def(i,1);
    i2 = local_edge_def(i,2);  
    for j = i:6
        j1 = local_edge_def(j,1);
        j2 = local_edge_def(j,2);

        Ae1(i,j) = Len(i)*Len(j)*(pxe*(c3(i1)*c4(i2)-c4(i1)*c3(i2))* ...
                                     (c3(j1)*c4(j2)-c4(j1)*c3(j2)) + ...
                                  pye*(c4(i1)*c2(i2)-c2(i1)*c4(i2))* ...
                                      (c4(j1)*c2(j2)-c2(j1)*c4(j2)) + ...
                                  pze*(c2(i1)*c3(i2)-c3(i1)*c2(i2))* ...
                                      (c2(j1)*c3(j2)-c3(j1)*c2(j2)));
        Ae1(j,i) = Ae1(i,j);

        % integral of the product of the scalar shape functions
        Ni1Nj1 = vole*NiNj(i1,j1);
        Ni1Nj2 = vole*NiNj(i1,j2);
        Ni2Nj1 = vole*NiNj(i2,j1);
        Ni2Nj2 = vole*NiNj(i2,j2);

        Ae2(i,j) = Len(i)*Len(j)*...
            (Ni1Nj1*(qxe*c2(i2)*c2(j2)+qye*c3(i2)*c3(j2)+qze*c4(i2)*c4(j2))- ...
             Ni1Nj2*(qxe*c2(i2)*c2(j1)+qye*c3(i2)*c3(j1)+qze*c4(i2)*c4(j1))- ...
             Ni2Nj1*(qxe*c2(i1)*c2(j2)+qye*c3(i1)*c3(j2)+qze*c4(i1)*c4(j2))+ ...
             Ni2Nj2*(qxe*c2(i1)*c2(j1)+qye*c3(i1)*c3(j1)+qze*c4(i1)*c4(j1)));  
                              
        Ae2(j,i) = Ae2(i,j);
 
    end
end
Ae1 = 4*Ae1*vole / (6*vole)^4;
Ae2 = Ae2 / (6*vole)^2;
Ae = Ae1 + Ae2;

%% 
% Calculate be: volume integral of F.Ni
if ~((Fxe==0) & (Fye==0) & (Fze==0))
    for i = 1:6
        i1 = local_edge_def(i,1);
        i2 = local_edge_def(i,2);
        
        be(i,1) = (Fxe*(c2(i2)-c2(i1)) + ...
                   Fye*(c3(i2)-c3(i1)) + ...
                   Fze*(c4(i2)-c4(i1)))*Len(i)/24;
    end
end


%%
function [local_edge_def] = fun_local_edge_def(nodes, edge2node_conne)

local_edge_def = zeros(6,2);
for i = 1:6
    n1 = edge2node_conne(i, 1);
    n2 = edge2node_conne(i, 2);
    ln1 = find(nodes == n1);
    ln2 = find(nodes == n2);
    local_edge_def(i,1) = ln1;
    local_edge_def(i,2) = ln2;
end
