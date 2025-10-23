function [Ex, Ey] = Efield(e, x, y, conn, u)
% E-field calculation (triangular elements)
% INPUT:
% e   : element ID
% x,y : x and y coordinates of all nodes in the mesh
% conn: connectivity matrix
% u   : potential values at all nodes in the mesh
% OUTPUT:
% Ex  : x-component of the electric field within the e-th element
% Ey  : y-component of the electric field within the e-th element
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

Ne = 3; % number of nodes in an element

% nodes
n1 = conn(e,1);
n2 = conn(e,2);
n3 = conn(e,3);

% potential values at the nodes
Ve = [u(n1), u(n2), u(n3)];

% calculate delxdelpsi, delxdelnu, delydelpsi, delydelnu
delxdelpsi = x(n2)-x(n1);
delxdelnu = x(n3)-x(n1);
delydelpsi = y(n2)-y(n1);
delydelnu = y(n3)-y(n1);

% calculate the determinant of the Jacobian matrix
Jdet = delxdelpsi*delydelnu-delxdelnu*delydelpsi;

% calculate the terms of inverse jacobian matrix
Jinv11 = delydelnu./Jdet;
Jinv12 = -delydelpsi./Jdet;
Jinv21 = -delxdelnu./Jdet;
Jinv22 = delxdelpsi./Jdet;

% calculate delNidelpsi, delNidelnu, delNjdelpsi, delNjdelnu
delNidelpsi = [-1 1 0];
delNidelnu  = [-1 0 1];

delVdelx = 0; delVdely = 0;
for i = 1:Ne
   delNidelx = Jinv11.*delNidelpsi(i)+Jinv12.*delNidelnu(i);
   delNidely = Jinv21.*delNidelpsi(i)+Jinv22.*delNidelnu(i);
   
   delVdelx = delVdelx + Ve(i)*delNidelx;
   delVdely = delVdely + Ve(i)*delNidely;
end  

Ex = -delVdelx;
Ey = -delVdely;

