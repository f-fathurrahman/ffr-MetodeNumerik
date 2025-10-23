function kernel = integrnd_quad_be(ksi, nu)
% Integrand for element right hand side vector (quadrilateral elements)
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
global i; 
global feq; global ee;

n1 = connq(ee,1); n2 = connq(ee,2); n3 = connq(ee,3); n4 = connq(ee,4);
x1 = xq(n1); x2 = xq(n2); x3 = xq(n3); x4 = xq(n4);
y1 = yq(n1); y2 = yq(n2); y3 = yq(n3); y4 = yq(n4);

ksii = [-1 1 1 -1];
nui  = [-1 -1 1 1];

% calculate delxdelksi, delxdelnu, delydelksi, delydelnu
T1 = -x1+x2+x3-x4;
T2 = x1-x2+x3-x4;
T3 = -x1-x2+x3+x4;
T4 = -y1+y2+y3-y4;
T5 = y1-y2+y3-y4;
T6 = -y1-y2+y3+y4;

delxdelksi = 0.25*(T1+nu*T2);
delxdelnu  = 0.25*(T3+ksi*T2);
delydelksi = 0.25*(T4+nu*T5);
delydelnu  = 0.25*(T6+ksi*T5);

% calculate the determinant of the Jacobian matrix
Jdet = delxdelksi.*delydelnu-delxdelnu.*delydelksi;

% calculate Ni
Ni = 0.25*(1+ksi.*ksii(i)).*(1+nu.*nui(i));

kernel = Ni.*feq.*Jdet;
