function kernel = integrnd_hexa_vole(ksi, nu, zeta)
% Integrand for the element area of hexahedral element
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
global ee;

n = zeros(1,8); x = zeros(1,8); y = zeros(1,8); z = zeros(1,8);
for ii = 1:8
    n(ii) = connq(ee,ii);
    x(ii) = xq(n(ii)); y(ii) = yq(n(ii)); z(ii) = zq(n(ii));
end

ksii = [-1 1 1 -1 -1 1 1 -1];
nui = [-1 -1 1 1 -1 -1 1 1];
zetai = [-1 -1 -1 -1 1 1 1 1];

delxdelksi = 0.125*sum(x.*ksii.*(1+nui.*nu).*(1+zetai.*zeta));
delxdelnu  = 0.125*sum(x.*nui.*(1+ksii.*ksi).*(1+zetai.*zeta));
delxdelzeta  = 0.125*sum(x.*zetai.*(1+ksii.*ksi).*(1+nui.*nu));
delydelksi = 0.125*sum(y.*ksii.*(1+nui.*nu).*(1+zetai.*zeta));
delydelnu  = 0.125*sum(y.*nui.*(1+ksii.*ksi).*(1+zetai.*zeta));
delydelzeta  = 0.125*sum(y.*zetai.*(1+ksii.*ksi).*(1+nui.*nu));
delzdelksi = 0.125*sum(z.*ksii.*(1+nui.*nu).*(1+zetai.*zeta));
delzdelnu  = 0.125*sum(z.*nui.*(1+ksii.*ksi).*(1+zetai.*zeta));
delzdelzeta  = 0.125*sum(z.*zetai.*(1+ksii.*ksi).*(1+nui.*nu));

J = [delxdelksi delydelksi delzdelksi; ...
    delxdelnu delydelnu delzdelnu; ...
    delxdelzeta delydelzeta delzdelzeta];

% calculate the determinant of the Jacobian matrix
Jdet = det(J);

kernel = Jdet;
