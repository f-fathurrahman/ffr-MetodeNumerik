function [anx, any] = unit_vector(x, y, node1, node2, nodeout)
% Calculate unit normal vector 
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

if (y(node1) < y(node2)+1e-8) && (y(node1) > y(node2)-1e-8) % ynode1=ynode2
   if y(nodeout) > y(node1)
      anx = 0; 
      any = 1;
   else 
      anx = 0; 
      any = -1;
   end
   
elseif (x(node1) < x(node2)+1e-8) && (x(node1) > x(node2)-1e-8) % xnode1=xnode2
   if x(nodeout) > x(node1)
      anx = 1; 
      any = 0;
   else 
      anx = -1; 
      any = 0;
   end
   
else
   theta = atan2(y(node2)-y(node1), x(node2)-x(node1));
   slope = tan(theta);
   ymid = 0.5*(y(node1)+y(node2));
   xmid = 0.5*(x(node1)+x(node2));

   term = -(x(node2)-x(node1)) / (y(node2)-y(node1));
   xu = (ymid-y(nodeout)+slope*x(nodeout)-xmid*term) / (slope-term);
   yu = y(nodeout) + slope*(xu - x(nodeout));

   thetap = atan2(yu-ymid, xu-xmid);
   anx = cos(thetap);
   any = sin(thetap);
end