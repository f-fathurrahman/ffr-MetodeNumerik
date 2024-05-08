%                       ***** PROGRAM 1.3 *****                            

% SIMULATION OF THE VERTICAL MOTION OF AN AIRFOIL IN WINDTUNNEL
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clear all;
format short e;

global ALPHA0 BETA PI U

% ***** ASSIGN VALUES TO ALPHA0, BETA, DT AND TMAX. SPECIFY INITIAL 
% POSITION AND VELOCITY OF THE AIRFOIL, THEN COMPUTE THE MOTION AT VARIOUS 
% WIND SPEEDS, HEADINGS ARE PRINTED FOR EACH WIND SPEED. *****

PI = 3.14159; 
% PI = pi;
ALPHA0 = 10.*PI/180.;
BETA = 0.00183;
DT = 0.1;
TMAX = 20.;
T0 = 0.0;
Z0 = 0.0;
V0 = 0.0;
U = 0.0;
DU = 100.0;



U = U + DU;

disp(['                 WIND SPEED U = ' num2str(U)]);

disp('    T              Z             V            THETA(DEG)   ');

%[T,SOL]=ode23(@F,[T0:DT:TMAX],[Z0 V0]);
[T,SOL]=ode45(@F,[T0:DT:TMAX],[Z0 V0]);
Z=SOL(:,1);
V=SOL(:,2);

ALPHA = ALPHA0 - atan(V/U);
ALPHAD = ALPHA*180./PI;

disp([T Z V ALPHAD]);

clear plot
plot(T,Z)
grid on
xlabel('T (nondimensional)'),ylabel('Z (nondimensional)')
title('Displacament of a wing with velocity U = 100 m/s')
