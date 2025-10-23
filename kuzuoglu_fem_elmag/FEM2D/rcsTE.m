function [rcs, farf] = rcsTE(Hz, conn, x, y, bound_node, bound_elm, scat_type)
% It calculates RCS for TE mode
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

global nu0;
global e0;
global k0;
global omg;
global rfar;

er = 1; % relative permittivity of the outer region (free-space)

rcs_term = zeros(1,360);
rcs_term2 = zeros(1,360);

snodes = zeros(1,length(x));
snodes(bound_node) = 1;
        
for i = 1:length(bound_elm)
    e = bound_elm(i);
    clear node1; clear node2; clear nodeout;
    
    n1 = conn(e,1);  n2 = conn(e,2);  n3 = conn(e,3);
            
    if (snodes(n1) && snodes(n2) && ~snodes(n3))
        node1 = n1;  node2 = n2;  nodeout = n3;
    elseif (~snodes(n1) && snodes(n2) && snodes(n3))
        node1 = n2;  node2 = n3;  nodeout = n1;
    elseif (snodes(n1) && ~snodes(n2) && snodes(n3))
        node1 = n3;  node2 = n1;  nodeout = n2;
    else
        disp('Error in the RCS calculation.');
        return;
    end
                
    delc = sqrt((x(node1)-x(node2))^2+(y(node1)-y(node2))^2);
    ym = 0.5*(y(node1)+y(node2));
    xm = 0.5*(x(node1)+x(node2));
    Hzm = 0.5*(Hz(node1)+Hz(node2));
        
    if clockwise([x(node1) x(node2) x(nodeout)], [y(node1) y(node2) y(nodeout)]) > 0
        [node2, node1] = deal(node1, node2); % swap node1 and node2
    end
    [anx, any] = unit_vector(x, y, node1, node2, nodeout);    
    
    Mz = (1/(1j*omg*e0*er))*normal_derivative(e, anx, any, x, y, conn, Hz);
    
    for phi = 0:359
        phir = phi*pi/180;
        Hf = exp(1j*k0*(xm*cos(phir)+ym*sin(phir)));
        
        Jxy = (any*sin(phir) + anx*cos(phir))*Hzm;
        
        rcs_term(phi+1) = rcs_term(phi+1) + Jxy*delc*Hf;
        
        if strcmpi(scat_type, 'diel')
            rcs_term2(phi+1) = rcs_term2(phi+1) + Mz*delc*Hf;
        end
    end
    
end

farf1 = rcs_term*sqrt(k0/(8*pi*rfar))*exp(-1j*k0*rfar+1j*pi/4);
farf2 = rcs_term2*(1/nu0)*sqrt(k0/(8*pi*rfar))*exp(-1j*k0*rfar-1j*3*pi/4);
farf = farf1+farf2;
rcs = 2*pi*rfar*abs(farf).^2;


%%***********************************************************************
%%***********************************************************************
%%***********************************************************************
function [area] = clockwise(svco_x, svco_y)
% if area is positive, triangle is oriented in clockwise direction,
% if negative, triangle is oriented in counter-clockwise direction

NV = length(svco_x);

area = 0;
for i = 1:NV
    if i == NV
        area = area + 0.5*(svco_x(1)-svco_x(i))*(svco_y(1)+svco_y(i));
    else
        area = area + 0.5*(svco_x(i+1)-svco_x(i))*(svco_y(i+1)+svco_y(i));
    end
end

%************************************************************************
function [ud] = normal_derivative(e, anx, any, x, y, conn, u)
% It calculates the normal derivative of a field using triangular elements

Ne = 3; % number of nodes in an element

% nodes
n1 = conn(e,1); n2 = conn(e,2); n3 = conn(e,3);

% field values at the nodes
uv = [u(n1), u(n2), u(n3)];

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
delNijdelpsi = [-1 1 0];
delNijdelnu  = [-1 0 1];

delEzdelx = 0;
delEzdely = 0;
for i = 1:Ne
   delNidelpsi = delNijdelpsi(i);
   delNidelnu  = delNijdelnu(i);

   delNidelx = Jinv11.*delNidelpsi+Jinv12.*delNidelnu;
   delNidely = Jinv21.*delNidelpsi+Jinv22.*delNidelnu;
   
   delEzdelx = delEzdelx + uv(i)*delNidelx;
   delEzdely = delEzdely + uv(i)*delNidely;
end  

ud = anx*delEzdelx + any*delEzdely;
