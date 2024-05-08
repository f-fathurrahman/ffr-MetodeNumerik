%                       ***** PROGRAM 1.6 *****                            

% GRAPHICAL PRESENTATION OF THE PATHS OF A GLIDER
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clf;
clear all;
format short e;

global A B


% ***** ASSIGN VALUES TO A, B, AND TO TIME INCREMENT DT *****

A = 1.5;
B = 0.06;
DT = 0.05;


% ***** COMPUTE TWO FLIGHT PATHS STARTING FROM THE ORIGIN WITH THE SAME
% INITIAL VELOCITY BUT IN DIFFERENT DIRECTIONS *****

W0 = 1.;
PI = 3.14159; 
% PI = pi;

for M=1:2

    if (M == 1) 
        THETA0(M) = -PI/2.;
    elseif (M == 2) 
        THETA0(M)=PI;
    end

    U0 = W0 * cos( THETA0(M) );
    V0 = W0 * sin( THETA0(M) );


% ***** INTEGRATE EQUATIONS (1.7.4) AND (1.7.5) USING THE FOURTH ORDER 
% RUNGE-KUTTA METHOD. X1, Y1 ARE RESPECTIVELY THE HORIZONTAL AND VERTICAL 
% COORDINATES OF THE GLIDER FOR M=1, AND X2, Y2 ARE THOSE FOR M=2 *****

    T = 0.;
    X = 0.;
    Y = 0.;
    U = U0;
    V = V0;
    for N=1:200 

        if ( M == 1 )
            X1(N) = X;
            Y1(N) = Y;
        elseif ( M == 2 )
            X2(N) = X;
            Y2(N) = Y;
        end

        [TSOL,SOL] = ode45(@F12, [T, T+DT], [X, U, Y, V]);
        T = TSOL(end);
        X = SOL(end,1);
        U = SOL(end,2);
        Y = SOL(end,3);
        V = SOL(end,4);

    end

end

plot(X1,Y1,'-',X2,Y2,'-.');
legend([['THETA0 = ' num2str(THETA0(1)*180./PI)] ' degrees'],[['THETA0 = ' num2str(THETA0(2)*180./PI)] ' degrees']);
title([[['Flight path of a glider starting from the origin having A = ' num2str(A)] [', B = ' num2str(B)]] [', W0 = ' num2str(W0)]] );

