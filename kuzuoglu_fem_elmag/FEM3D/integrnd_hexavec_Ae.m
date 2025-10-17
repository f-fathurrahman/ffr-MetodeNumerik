function kernel = integrnd_hexavec_Ae(ksi, nu, zeta)
% Integrand for the element matrix of hexahedral edge element
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
global elm2edge_connq; global Lenq;

if (i >= 1) & (i <= 4)
    iax = 1; % parallel to ksi axis
elseif (i >= 5) & (i <= 8)
    iax = 2; % parallel to nu axis
elseif (i >= 9) & (i <= 12)
    iax = 3; % parallel to zeta axis
end
    
if (j >= 1) & (j <= 4)
    jax = 1; % parallel to ksi axis
elseif (j >= 5) & (j <= 8)
    jax = 2; % parallel to nu axis
elseif (j >= 9) & (j <= 12)
    jax = 3; % parallel to zeta axis
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

T1 = Jinv(2,1).*Jinv(3,2) - Jinv(3,1).*Jinv(2,2);
T2 = Jinv(3,1).*Jinv(1,2) - Jinv(1,1).*Jinv(3,2);
T3 = Jinv(1,1).*Jinv(2,2) - Jinv(2,1).*Jinv(1,2);
T4 = Jinv(2,1).*Jinv(3,3) - Jinv(3,1).*Jinv(2,3);
T5 = Jinv(3,1).*Jinv(1,3) - Jinv(1,1).*Jinv(3,3);
T6 = Jinv(1,1).*Jinv(2,3) - Jinv(2,1).*Jinv(1,3);
T7 = Jinv(2,2).*Jinv(3,3) - Jinv(3,2).*Jinv(2,3);
T8 = Jinv(3,2).*Jinv(1,3) - Jinv(1,2).*Jinv(3,3);
T9 = Jinv(1,2).*Jinv(2,3) - Jinv(2,2).*Jinv(1,3);

eksii = [0 0 0 0 -1 -1 1 1 -1 1 -1 1];
enui = [-1 1 -1 1 0 0 0 0 -1 -1 1 1];
ezetai = [-1 -1 1 1 -1 1 -1 1 0 0 0 0];

R1i = (Len(i)/8) .* (1+nu.*enui(i)) .* ezetai(i);
R2i = (Len(i)/8) .* (1+zeta.*ezetai(i)) .* enui(i);
R3i = (Len(i)/8) .* (1+ksi.*eksii(i)) .* ezetai(i);
R4i = (Len(i)/8) .* (1+zeta.*ezetai(i)) .* eksii(i);
R5i = (Len(i)/8) .* (1+ksi.*eksii(i)) .* enui(i);
R6i = (Len(i)/8) .* (1+nu.*enui(i)) .* eksii(i);

R1j = (Len(j)/8) .* (1+nu.*enui(j)) .* ezetai(j);
R2j = (Len(j)/8) .* (1+zeta.*ezetai(j)) .* enui(j);
R3j = (Len(j)/8) .* (1+ksi.*eksii(j)) .* ezetai(j);
R4j = (Len(j)/8) .* (1+zeta.*ezetai(j)) .* eksii(j);
R5j = (Len(j)/8) .* (1+ksi.*eksii(j)) .* enui(j);
R6j = (Len(j)/8) .* (1+nu.*enui(j)) .* eksii(j);

if (iax == 1) & (jax == 1) % ksi-ksi
    kernel1 = pxeq.*(R1i.*T4+R2i.*T1).*(R1j.*T4+R2j.*T1) + ...
              pyeq.*(R1i.*T5+R2i.*T2).*(R1j.*T5+R2j.*T2) + ...
              pzeq.*(R1i.*T6+R2i.*T3).*(R1j.*T6+R2j.*T3);
    kernel2 = (Len(i).*Len(j)/64).*(1+nu.*enui(i)).*(1+zeta.*ezetai(i)).* ...
        (1+nu.*enui(j)).*(1+zeta.*ezetai(j)) .* ...
        (qxeq.*Jinv(1,1).^2+qyeq.*Jinv(2,1).^2+qzeq.*Jinv(3,1).^2);
    
elseif (iax == 2) & (jax == 2) % nu-nu
    kernel1 = pxeq.*(R4i.*T1-R3i.*T7).*(R4j.*T1-R3j.*T7) + ...
              pyeq.*(R4i.*T2-R3i.*T8).*(R4j.*T2-R3j.*T8)+ ...
              pzeq.*(R4i.*T3-R3i.*T9).*(R4j.*T3-R3j.*T9);
    kernel2 = (Len(i).*Len(j)/64).*(1+ksi.*eksii(i)).*(1+zeta.*ezetai(i)).* ...
        (1+ksi.*eksii(j)).*(1+zeta.*ezetai(j)) .* ...
        (qxeq.*Jinv(1,2).^2+qyeq.*Jinv(2,2).^2+qzeq.*Jinv(3,2).^2);
    
elseif (iax == 3) & (jax == 3) % zeta-zeta
    kernel1 = pxeq.*(R5i.*T7+R6i.*T4).*(R5j.*T7+R6j.*T4) + ...
              pyeq.*(R5i.*T8+R6i.*T5).*(R5j.*T8+R6j.*T5)+ ...
              pzeq.*(R5i.*T9+R6i.*T6).*(R5j.*T9+R6j.*T6);
    kernel2 = (Len(i).*Len(j)/64).*(1+ksi.*eksii(i)).*(1+nu.*enui(i)).* ...
        (1+ksi.*eksii(j)).*(1+nu.*enui(j)) .* ...
        (qxeq.*Jinv(1,3).^2+qyeq.*Jinv(2,3).^2+qzeq.*Jinv(3,3).^2);
    
elseif (iax == 1) & (jax == 2) % ksi-nu
    kernel1 = -pxeq.*(R1i.*T4+R2i.*T1).*(R4j.*T1-R3j.*T7) - ...
               pyeq.*(R1i.*T5+R2i.*T2).*(R4j.*T2-R3j.*T8) - ...
               pzeq.*(R1i.*T6+R2i.*T3).*(R4j.*T3-R3j.*T9);
    kernel2 = (Len(i).*Len(j)/64).*(1+nu.*enui(i)).*(1+zeta.*ezetai(i)).* ...
        (1+ksi.*eksii(j)).*(1+zeta.*ezetai(j)) .* ...
        (qxeq.*Jinv(1,1).*Jinv(1,2)+qyeq.*Jinv(2,1).*Jinv(2,2)+ ...
        qzeq.*Jinv(3,1).*Jinv(3,2));
    
elseif (iax == 2) & (jax == 1) % nu-ksi
    kernel1 = -pxeq.*(R4i.*T1-R3i.*T7).*(R1j.*T4+R2j.*T1) - ...
               pyeq.*(R4i.*T2-R3i.*T8).*(R1j.*T5+R2j.*T2)- ...
               pzeq.*(R4i.*T3-R3i.*T9).*(R1j.*T6+R2j.*T3);
    kernel2 = (Len(i).*Len(j)/64).*(1+nu.*enui(j)).*(1+zeta.*ezetai(j)).* ...
        (1+ksi.*eksii(i)).*(1+zeta.*ezetai(i)) .* ...
        (qxeq.*Jinv(1,1).*Jinv(1,2)+qyeq.*Jinv(2,1).*Jinv(2,2)+ ...
        qzeq.*Jinv(3,1).*Jinv(3,2));
    
elseif (iax == 1) & (jax == 3) % ksi-zeta
    kernel1 = -pxeq.*(R1i.*T4+R2i.*T1).*(R5j.*T7+R6j.*T4) - ...
               pyeq.*(R1i.*T5+R2i.*T2).*(R5j.*T8+R6j.*T5) - ...
               pzeq.*(R1i.*T6+R2i.*T3).*(R5j.*T9+R6j.*T6);
    kernel2 = (Len(i).*Len(j)/64).*(1+nu.*enui(i)).*(1+zeta.*ezetai(i)).* ...
        (1+ksi.*eksii(j)).*(1+nu.*enui(j)) .* ...
        (qxeq.*Jinv(1,1).*Jinv(1,3)+qyeq.*Jinv(2,1).*Jinv(2,3)+ ...
        qzeq.*Jinv(3,1).*Jinv(3,3));
    
elseif (iax == 3) & (jax == 1) % zeta-ksi
    kernel1 = -pxeq.*(R5i.*T7+R6i.*T4).*(R1j.*T4+R2j.*T1) - ...
               pyeq.*(R5i.*T8+R6i.*T5).*(R1j.*T5+R2j.*T2) - ...
               pzeq.*(R5i.*T9+R6i.*T6).*(R1j.*T6+R2j.*T3);
    kernel2 = (Len(i).*Len(j)/64).*(1+nu.*enui(j)).*(1+zeta.*ezetai(j)).* ...
        (1+ksi.*eksii(i)).*(1+nu.*enui(i)) .* ...
        (qxeq.*Jinv(1,1).*Jinv(1,3)+qyeq.*Jinv(2,1).*Jinv(2,3)+ ...
        qzeq.*Jinv(3,1).*Jinv(3,3));
    
elseif (iax == 2) & (jax == 3) % nu-zeta
    kernel1 = pxeq.*(R4i.*T1-R3i.*T7).*(R5j.*T7+R6j.*T4) + ...
              pyeq.*(R4i.*T2-R3i.*T8).*(R5j.*T8+R6j.*T5)+ ...
              pzeq.*(R4i.*T3-R3i.*T9).*(R5j.*T9+R6j.*T6);
    kernel2 = (Len(i).*Len(j)/64).*(1+ksi.*eksii(i)).*(1+zeta.*ezetai(i)).* ...
        (1+ksi.*eksii(j)).*(1+nu.*enui(j)) .* ...
        (qxeq.*Jinv(1,2).*Jinv(1,3)+qyeq.*Jinv(2,2).*Jinv(2,3)+ ...
        qzeq.*Jinv(3,2).*Jinv(3,3));
    
elseif (iax == 3) & (jax == 2) % zeta-nu
    kernel1 = pxeq.*(R5i.*T7+R6i.*T4).*(R4j.*T1-R3j.*T7) + ...
              pyeq.*(R5i.*T8+R6i.*T5).*(R4j.*T2-R3j.*T8)+ ...
              pzeq.*(R5i.*T9+R6i.*T6).*(R4j.*T3-R3j.*T9);
    kernel2 = (Len(i).*Len(j)/64).*(1+ksi.*eksii(j)).*(1+zeta.*ezetai(j)).* ...
        (1+ksi.*eksii(i)).*(1+nu.*enui(i)) .* ...
        (qxeq.*Jinv(1,2).*Jinv(1,3)+qyeq.*Jinv(2,2).*Jinv(2,3)+ ...
        qzeq.*Jinv(3,2).*Jinv(3,3));
end

kernel = (kernel1 + kernel2).*Jdet;
     