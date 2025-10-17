function kernel = integrnd_hexavec_be(ksi, nu, zeta)
% Integrand for the element right hand side vector of hexahedral edge element
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
global i; global ee;
global Fxeq; global Fyeq; global Fzeq;
global elm2edge_connq; global Lenq;

if (i >= 1) & (i <= 4)
    iax = 1; % parallel to ksi axis
elseif (i >= 5) & (i <= 8)
    iax = 2; % parallel to nu axis
elseif (i >= 9) & (i <= 12)
    iax = 3; % parallel to zeta axis
end

edges = elm2edge_connq(ee,:);
Len = Lenq(edges);
n = connq(ee,:);
x = xq(n); y = yq(n); z = zq(n); 

ksii = [-1 1 1 -1 -1 1 1 -1];
nui = [-1 -1 1 1 -1 -1 1 1];
zetai = [-1 -1 -1 -1 1 1 1 1];

delxdelksi = 0.125*sum(x.*ksii.*(1+nui.*nu).*(1+zetai.*zeta));
delxdelnu = 0.125*sum(x.*nui.*(1+ksii.*ksi).*(1+zetai.*zeta));
delxdelzeta = 0.125*sum(x.*zetai.*(1+ksii.*ksi).*(1+nui.*nu));
delydelksi = 0.125*sum(y.*ksii.*(1+nui.*nu).*(1+zetai.*zeta));
delydelnu = 0.125*sum(y.*nui.*(1+ksii.*ksi).*(1+zetai.*zeta));
delydelzeta = 0.125*sum(y.*zetai.*(1+ksii.*ksi).*(1+nui.*nu));
delzdelksi = 0.125*sum(z.*ksii.*(1+nui.*nu).*(1+zetai.*zeta));
delzdelnu = 0.125*sum(z.*nui.*(1+ksii.*ksi).*(1+zetai.*zeta));
delzdelzeta = 0.125*sum(z.*zetai.*(1+ksii.*ksi).*(1+nui.*nu));

J = [delxdelksi delydelksi delzdelksi; ...
    delxdelnu delydelnu delzdelnu; ...
    delxdelzeta delydelzeta delzdelzeta];

% calculate the determinant of the Jacobian matrix
Jdet = det(J);

% calculate the terms of inverse jacobian matrix
Jinv = inv(J);

eksii = [0 0 0 0 -1 -1 1 1 -1 1 -1 1];
enui = [-1 1 -1 1 0 0 0 0 -1 -1 1 1];
ezetai = [-1 -1 1 1 -1 1 -1 1 0 0 0 0];

if (iax == 1) % parallel to ksi axis
    kernel = (Len(i)/8) .* (1+nu.*enui(i)) .* (1+zeta.*ezetai(i)) .* ...
             (Jinv(1,1).*Fxeq + Jinv(2,1).*Fyeq + Jinv(3,1).*Fzeq);
elseif (iax == 2) % parallel to nu axis
    kernel = (Len(i)/8) .* (1+ksi.*eksii(i)) .* (1+zeta.*ezetai(i)) .* ...
             (Jinv(1,2).*Fxeq + Jinv(2,2).*Fyeq + Jinv(3,2).*Fzeq);
elseif (iax == 3) % parallel to zeta axis
    kernel = (Len(i)/8) .* (1+ksi.*eksii(i)) .* (1+nu.*enui(i)) .* ...
             (Jinv(1,3).*Fxeq + Jinv(2,3).*Fyeq + Jinv(3,3).*Fzeq);
end

kernel = kernel.*Jdet;