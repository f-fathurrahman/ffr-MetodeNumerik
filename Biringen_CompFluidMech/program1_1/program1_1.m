%                       ***** PROGRAM 1.1 *****                            

% COMPARISON OF THE MOTION OF A STEEL SPHERE  DROPPED IN AIR WITH THAT IN
% VACUUM, BASED ON FOURTH-ORDER RUNGE-KUTTA METHOD
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clear all; close all;
%format short e;

% Global variables !!!
global A B C D NU

% Data
T0 = 0.;
Z0 = 0.;
V0 = 0.;
RHO = 8000.0;
RHOF = 1.22;
G = 9.8;
H = 0.1;
TMAX = 10.0;
RHOBAR = RHOF/RHO;
NU = 0.0000149;
A = 1.+RHOBAR/2.;
B = (1.-RHOBAR)*G;
D = 0.01;
C = 3.*RHOBAR/(4.*D);

% Initialization
% FIXME: Arrays should be pre-allocated here
T = T0;
Z = Z0;
V = V0;
RE = V*D/NU; % Reynold's number

% Header
disp('    T          Z(IN VACUUM)      Z       V(IN VACUUM)      V            RE');
disp('  (SEC)             (M)         (M)         (M/SEC)     (M/SEC)           ');
N = 1;

while (T<=TMAX)

    % The position and velocity in vacuum (using analytic expression)
    ZV = Z0 + V0*T + G*T*T/2.0;
    VV = V0 + G*T;

    disp([T, ZV, Z, VV, V, RE]);

    % Runge-Kutta method begins here
    % 
    D1Z = H*V;
    D1V = H*F(V);
    %
    D2Z = H*(V + D1V/2.);
    D2V = H*F(V + D1V/2.);
    %
    D3Z = H*(V + D2V/2.);
    D3V = H*F(V + D2V/2.);
    %
    D4Z = H*(V + D3V);
    D4V = H*F(V + D3V);
    %
    T = T + H; % increment time
    N = N + 1; % array index: FIXME: remove this, use preallocated array
    #
    Z = Z + (D1Z + 2.*D2Z + 2.*D3Z + D4Z)/6.;
    V = V + (D1V + 2.*D2V + 2.*D3V + D4V)/6.;
    RE = V*D/NU;
    %
    % FIXME: These arrays will grow in size
    %
    TT(N) = T;
    ZZ(N) = Z;
    VV1(N) = V;
    %
    % Plotting stuffs
    %
    %clear plot;
    %plot(TT,ZZ);
    %grid on
    %xlabel('t (s)'),ylabel('z (m), v (m/s)');
    %title('Sphere diameter, D = 0.01');
    %legend(', z(t); -----, v(t)')
    %hold on
    %clear plot;
    %plot(TT, VV1,'--');

% ***** REPEAT THE COMPUTATION IF T<=TMAX, OTHERWISE STOP *****

end
