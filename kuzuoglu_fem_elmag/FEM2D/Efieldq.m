function [Ex, Ey] = Efieldq(e, x, y, conn, u)
% E-field calculation (quadrilateral elements)
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

Ne = 4; % number of nodes in an element

% nodes
n1 = conn(e,1); n2 = conn(e,2); n3 = conn(e,3); n4 = conn(e,4);
x1 = x(n1); x2 = x(n2); x3 = x(n3); x4 = x(n4);
y1 = y(n1); y2 = y(n2); y3 = y(n3); y4 = y(n4);

% potential values at the nodes
Ve = [u(n1), u(n2), u(n3), u(n4)];

nu = 0; psi = 0; % center of quad (derivative is evaluated at the midpoint)

psii = [-1 1 1 -1];
nui  = [-1 -1 1 1];

% calculate delxdelpsi, delxdelnu, delydelpsi, delydelnu
T1 = -x1+x2+x3-x4;
T2 = x1-x2+x3-x4;
T3 = -x1-x2+x3+x4;
T4 = -y1+y2+y3-y4;
T5 = y1-y2+y3-y4;
T6 = -y1-y2+y3+y4;

delxdelpsi = 0.25*(T1+nu*T2);
delxdelnu  = 0.25*(T3+psi*T2);
delydelpsi = 0.25*(T4+nu*T5);
delydelnu  = 0.25*(T6+psi*T5);

% calculate the determinant of the Jacobian matrix
Jdet = delxdelpsi.*delydelnu-delxdelnu.*delydelpsi;

% calculate the terms of inverse jacobian matrix
Jinv11 = delydelnu./Jdet;
Jinv12 = -delydelpsi./Jdet;
Jinv21 = -delxdelnu./Jdet;
Jinv22 = delxdelpsi./Jdet;

delVdelx = 0; delVdely = 0;
for i = 1:Ne
    delNidelpsi = 0.25*psii(i)*(1+nu.*nui(i));
    delNidelnu  = 0.25*nui(i)*(1+psi.*psii(i));

    delNidelx = Jinv11.*delNidelpsi+Jinv12.*delNidelnu;
    delNidely = Jinv21.*delNidelpsi+Jinv22.*delNidelnu;
   
    delVdelx = delVdelx + Ve(i)*delNidelx;
    delVdely = delVdely + Ve(i)*delNidely;
end    

Ex = -delVdelx;
Ey = -delVdely;

