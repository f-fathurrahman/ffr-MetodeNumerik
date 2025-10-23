function [Ae, be, areae] = element_matrix_tri2(e, x, y, conn, pxe, pye, qe, fe)
% Element matrix equation with numerical calculation (triangular elements)
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

Ne = 3; % number of nodes in an element
Ae1 = zeros(Ne,Ne);
Ae2 = zeros(Ne,Ne);
be = zeros(Ne,1);

n1 = conn(e,1);  n2 = conn(e,2);  n3 = conn(e,3); % nodes of the element
x1 = x(n1);  x2 = x(n2);  x3 = x(n3); % x-coordinates of the nodes
y1 = y(n1);  y2 = y(n2);  y3 = y(n3); % y-coordinates of the nodes

% Jacobian entries
delxdelpsi = x2-x1;
delxdelnu = x3-x1;
delydelpsi = y2-y1;
delydelnu = y3-y1;

Jdet = delxdelpsi*delydelnu-delxdelnu*delydelpsi; %determinant of jacobian
areae = 0.5*Jdet; % area of the element

% entries of the inverse Jacobian matrix
Jinv11 = delydelnu./Jdet;
Jinv12 = -delydelpsi./Jdet;
Jinv21 = -delxdelnu./Jdet;
Jinv22 = delxdelpsi./Jdet;

% calculate delNidelpsi, delNidelnu, delNjdelpsi, delNjdelnu
delNijdelpsi = [-1 1 0];
delNijdelnu  = [-1 0 1];

% Element matrix (first part)
for i = 1:Ne
    for j = i:Ne
        Ae1(i,j)=(pxe*(Jinv11.*delNijdelpsi(i)+Jinv12.*delNijdelnu(i)).* ...
                       (Jinv11.*delNijdelpsi(j)+Jinv12.*delNijdelnu(j))+ ...
                  pye*(Jinv21.*delNijdelpsi(i)+Jinv22.*delNijdelnu(i)).* ...
                       (Jinv21.*delNijdelpsi(j)+Jinv22.*delNijdelnu(j)))*Jdet*0.5;
        Ae1(j,i) = Ae1(i,j);
    end
end

% Element matrix (second part)

% Gauss weights and points (3-point)
gw = [1/6 1/6 1/6];   % Gauss weights 
ksi = [0.5 0 0.5];    % Gauss point ksi
nu = [0 0.5 0.5];     % Gauss point nu

for i = 1:Ne
    for j = i:Ne
        sumg = 0;
        for k = 1:length(gw)
            if i == 1
                Ni = 1-ksi(k)-nu(k);
            elseif i == 2
                Ni = ksi(k);
            else
                Ni = nu(k);
            end
            
            if j == 1
                Nj = 1-ksi(k)-nu(k);
            elseif j == 2
                Nj = ksi(k);
            else
                Nj = nu(k);
            end        
            
            sumg = sumg+gw(k)*qe*Ni*Nj*Jdet;
        end
        
        Ae2(i,j) = sumg;
        Ae2(j,i) = Ae2(i,j);
    end
end

% Element matrix
Ae = Ae1 + Ae2;

% Element right-hand-side vector
for i = 1:Ne
    sumg = 0;
    for k = 1:length(gw)
        if i == 1
            Ni = 1-ksi(k)-nu(k);
        elseif i == 2
            Ni = ksi(k);
        else
            Ni = nu(k);
        end
        
        sumg = sumg+gw(k)*fe*Ni*Jdet;
    end
    
    be(i) = sumg;
end
