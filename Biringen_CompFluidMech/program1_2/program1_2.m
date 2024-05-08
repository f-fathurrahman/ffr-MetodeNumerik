%                       ***** PROGRAM 1.2 *****                            

% MOTION OF A SIMPLE PENDULUM IN AIR
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clear all;
format short e;

global A B C D G L NU

% ***** SPECIFY DATA *****

RHO = 2500.;
RHOF = 1.22;
%PI = 3.14159; 
PI = pi;
DT = 0.1;
TMAX = 20.;
D = 0.01;
G = 9.8;
L = 2.0;
NU = 1.49E-5;
RHOBAR = RHOF/RHO;
A = 1.+RHOBAR/2.;
B = (1.-RHOBAR)*G;
C = 3.*RHOBAR/(4.*D);


% ***** INITIALIZE THE PROBLEM *****

T0 = 0.0;
V0 = 0.0;
THETA0 = 45.0;
Y0 = L*THETA0*PI/180.0;


% ***** INTEGRATE THE EQUATIONS OF MOTION BY USING FOURTH-ORDER RUNGE-KUTTA
% METHOD, STOP THE COMPUTATION WHEN T>TMAX *****


disp('                     IN  AIR      ');
disp('    T                 THETA     V      ');
disp('  (SEC)                 (DEG)   M/SEC)   ');

%[T,SOLN]=ode23(@FN,[T0:DT:TMAX],[Y0 V0]);
[T,SOLN]=ode45(@FN,[T0:DT:TMAX],[Y0 V0]);
Y=SOLN(:,1);
V=SOLN(:,2);
THETA = Y/L * 180.0/PI;

disp([T  THETA V]);

%PLOT THE DATA
clear plot;
plot(T,THETA)
grid on
xlabel('t (s)'),ylabel('\theta (deg)')
title('Angular displacement  for glass sphere in air')




