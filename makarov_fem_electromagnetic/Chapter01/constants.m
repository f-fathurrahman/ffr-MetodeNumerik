%   Compute electromagnetic constants
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%%  EM constants
eps0        = 8.85418782e-012;                  %   dielectric permittivity of vacuum(~air), F/m
mu0         = 1.25663706e-006;                  %   magnetic permeability of vacuum(~air), H/m
c0          = 1/sqrt(eps0*mu0);                 %   speed of light in vacuum(~air), m/s
eta0        = sqrt(mu0/eps0);                   %   vacuum/air impedance, ohms ~377 ohms