%                       ***** PROGRAM 1.5 *****                            

% FINDING THE MAXIMUM RANGE OF A STEEL SPHERICAL PROJECTILE
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clear all;
format short e;

global A B C D NU PI UF VF


% ***** ASSIGN INPUT DATA AND PRINT HEADINGS ON A NEW PAGE *****

RHO = 8000.0;
RHOF = 1.22;
G = 9.8;
DT = 0.1;
EPSILON = 0.001;
T0 = 0.0;
X0 = 0.0;
Y0 = 0.0;
W0 = 50.0;
RHOBAR = RHOF/RHO;
NU = 1.49E-5;
PI = 3.14159; 
% PI = pi;
disp([[[['CONSIDER A STEEL SPHERICAL PROJECTILE OF ' num2str(D)] ' M DIAMETER SHOOTING IN AIR WITH AN INITIAL SPEED OF '] num2str(W0)] 'M/SEC.']);
disp('SYMBOLS USED IN THE FOLLOWING TABLE ARE ;');
disp('UF = HORIZONTAL WIND SPEED');
disp('THETA = THE OPTIMUM SHOOTING ANGLE');
disp('XRMAX = THE MAXIMUM RANGE');
disp('N = NUMBER OF STEPS USED IN HALF-INTERNAL METHOD TO OBTAIN THE PRINTED RESULT');
disp([['THE ERROR IN THETA IS OF THE ORDER OF ' num2str(EPSILON)] 'DEGREES.']);

disp('      UF          THETA        XRMAX          N ');
disp('    (M/SEC)       (DEG)         (M)             ');
disp('    -------       -----        -----         ---')
D = 2.5E-3;

for J=1:41

A = 1. + RHOBAR/2.;
B = (1.-RHOBAR)*G;
C = 3.*RHOBAR/(4.*D);



% ***** ASSIGN  WIND CONDITIONS - NO WIND *****

VF = 0.;
UF = 0.;


    
% ***** START WITH A HORIZONTAL SHOOTING POSITION AND WITH A 10-DEGREE 
% INCREMENT *****

    THEOLD = 0.;
    XROLD = 40.;
    DELTA = 10.;
    N = 0;


% ***** ACCORDING TO HALF-INTERNAL METHOD, THE ABSCISSA OF THE NEXT POINT
% ON THE SHOOTING ANGLE-RANGE CURVE IS *****

    while (abs(DELTA) >= EPSILON)
        THENEW = THEOLD + DELTA;
        N = N + 1;


% ***** TO FIND THE ORDINATE OF THAT POINT, WE FIRST USE RUNGE-KUTTA METHOD
% TO DETERMINE THE POINTS P AND Q ON THE TRAJECTORY, THEN USE EQUATION
% (1.5.3) TO COMPUTE THE RANGE *****

        T = T0;
        X = X0;
        Y = Y0;
        U = W0 * cos(THENEW*PI/180.);
        V = W0 * sin(THENEW*PI/180.);

        while (Y >= 0.)
            [TSOL,SOL] = ode45(@FXY, [T, T+DT], [X, U, Y, V]);
            T = TSOL(end);
            X = SOL(end,1);
            U = SOL(end,2);
            Y = SOL(end,3);
            V = SOL(end,4);
        end
        XQ = X;
        YQ = Y;
        [TSOL,SOL] = ode45(@FXY, [T, T-DT], [X, U, Y, V]);
        T = TSOL(end);
        X = SOL(end,1);
        U = SOL(end,2);
        Y = SOL(end,3);
        V = SOL(end,4);
        XP = X;
        YP = Y;
        XRNEW = (XQ*YP-XP*YQ)/(YP-YQ);


% ***** CHANGE THE VALUE OF DELTA IF XROLD>=XRNEW. STOP HALF-INTERVAL 
% METHOD IF DELTA<EPSILON *****

        if (XROLD>=XRNEW)
            DELTA = -DELTA/2.;
        else
% ***** CALL THE NEW POINT AN OLD ONE, THEN LOCATE THE NEXT NEW POINT *****
            THEOLD = THENEW;
            XROLD = XRNEW;
        end
    end


% ***** THE MAXIMUM RANGE, XRMAX, AND THE OPTIMUM SHOOTING ANGLE, THETA,
% ARE DETERMINED *****

    if (XROLD>=XRNEW)
        XRMAX = XROLD;
        THETA = THEOLD;
    else
        XRMAX = XRNEW;
        THETA = THENEW;
    end

    disp([UF THETA XRMAX N]);
    VAR1(J)=D;
    VAR2(J)=THETA;   
        if D<=1.0E-2
                D=D+0.0005;  
            elseif D<=2.0E-2&D>1.0E-2
                D=D+0.00125;
            elseif D>2.0E-2   
                D=D+0.005;
        end
   
end
%PLOT RESULTS	

semilogx(VAR1,VAR2,'o')
grid on
xlabel('Diameter (m)'), ylabel(texlabel('theta (deg)'))
title('Optimum shooting angle, steel spherical projectile, still air, initial velocity 50 m/s')
