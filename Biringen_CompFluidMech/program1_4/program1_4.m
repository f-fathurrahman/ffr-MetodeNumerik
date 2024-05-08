%                       ***** PROGRAM 1.4 *****                            

% MOTION AND TRAJECTORY OF A SPHERICAL PROJECTILE IN A FLUID
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)
clc;
clear all;
format short e;

global A B C D NU PI

% ***** ASSIGN INPUT DATA AND START THE OUTPUT ON A NEW PAGE *****
RHO = 8000.0;
RHOF = 1.22;
G = 9.8;
DT = 0.1;
TMAX = 10;
IFRPRT = 5;
T0 = 0.0;
X0 = 0.0;
Y0 = 0.0;
W0 = 50.0;
RHOBAR = RHOF/RHO;
NU = 1.49E-5;
PI = 3.14159; 
% PI = pi;
D = 0.05;
A = 1. + RHOBAR/2.;
B = (1. - RHOBAR)*G;
C = 3.*RHOBAR/(4.*D);

% ***** PRINTING HEADINS AND ASSIGN INITIAL CONDITIONS FOR EACH OF THE 3
% CASES WITH THETA0 = 30, 45 AND 60 DEGREES, RESPECTIVELY. THEN PRINT THE
% INITIAL VALUES *****

THETA0 = 15.0;

    J=1
    THETA0 = THETA0 + 15.0;
    I = 1;
    T(I) = T0;
    X(I) = X0;
    Y(I) = Y0;
    U(I) = W0 * cos(THETA0*PI/180.);
    V(I) = W0 * sin(THETA0*PI/180.);

% ***** FOURTH-ORDER RUNGE-KUTTA METHOD IS USED TO INTEGRATE THE EQUATIONS
% OF MOTION, STOP THE COMPUTATION FOR ONE CASE IF THE PROJECTILE DROPS BACK
% TO ITS INITIAL LEVEL OR BELOW, PRINT DATA WHEN THE LAST TIME STEP IS
% REACHED OR WHEN THE COUNT I IS AN INTEGRAL MULTIPLE OF IFRPRT *****

    disp([[[[[['CASE(' num2str(J)] ') ***** INITIAL SPEED ='] num2str(W0)] ' M/SEC.    THETA0 ='] num2str(THETA0)] ' DEGREES']);
    disp('       T            X            Y            U            V   ');
    disp('     (SEC)         (M)          (M)        (M/SEC)      (M/SEC)');

    disp([T(I) X(I) Y(I) U(I) V(I)]);

    while (Y(I) >= 0.)
        [TSOL,SOL] = ode45(@FXY, [T(I), T(I)+DT], [X(I), U(I), Y(I), V(I)]);
        I = I + 1;
        T(I)=TSOL(end);
        X(I)=SOL(end,1);
        U(I)=SOL(end,2);
        Y(I)=SOL(end,3);
        V(I)=SOL(end,4);
        if (rem(I-1,IFRPRT) == 0)
            disp([T(I) X(I) Y(I) U(I) V(I)]);
        end
    end
    disp([T(I) X(I) Y(I) U(I) V(I)]);
    clear plot
    plot(X,Y)
    xlabel('X (m)'),ylabel('Y (m)')
    title('Trajectory of a steel sphere in air, initial velocity 50 m/s; initial projectile elevation 30 degrees')

