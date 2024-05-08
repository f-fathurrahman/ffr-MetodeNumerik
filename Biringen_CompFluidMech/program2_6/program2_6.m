%                       ***** PROGRAM 2.6 *****                            

% NUMERICAL SOLUTION OF FLOW THROUGH A CHAMBER HAVING AN IRREGULAR BOUNDARY
% USING SUCCESSIVE OVERRELAXATION METHOD
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clf;
clear all;
format short e;

global M N H JLAST;

M = 25;
N = 17;

SYMBOL = zeros(10,1);
CONST = zeros(10,1);
JLAST = zeros(M,1);
PSI = zeros(M,N);
X = zeros(M,1);
Y = zeros(N,1);


% ***** ASSIGN INPUT DATA *****

H = 0.25;
SYMBOL = '+osdv^x*ph.';
CONST = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0];


% ***** DEFINE THE LAST VALUE OF J ALONG EACH VERTICAL GRID LINE *****

for I=1:N
    JLAST(I) = N;
end

for I=18:M
    JLAST(I) = N - I + N;
end


% ***** ASSIGN BOUNDARY VALUES TO PSI *****

for J=1:N
    PSI(1,J) = 1.0;
end

for J=1:9
    PSI(M,J) = 1.0;
end
   
for I=2:M-1
    if (I > 7 && I < 22)
        PSI(I,1) = 0.0;
    else
        PSI(I,1) = 1.0;
    end
    J = JLAST(I);
    PSI(I,J) = 1.0;
end


% ***** ASSIGN GUESSED VALUES TO PSI *****

for I=2:M-1
    JEND = JLAST(I) - 1;
    for J=2:JEND
        PSI(I,J) = 0.5;
    end
end


% ***** ITERATE THE SUCCESSIVE OVERRELAXATION FORMULA UNTIL THE SUM OF
% ERROS IS LESS THAN OR EQUAL TO 0.001 *****

ALPHA = cos(pi/M) + cos(pi/N);
OMEGA = (8.0-4.0*sqrt(4.0-ALPHA^2)) / ALPHA^2;
ITER = 0;
ERROR = 1.0;
while (ERROR > 0.001)
    ITER = ITER + 1;
    ERROR = 0.0;
    for I=2:M-1
        JEND = JLAST(I) - 1;
        for J=2:JEND
            PSIPRV = PSI(I,J);
            PSI(I,J) = (1.0-OMEGA)*PSI(I,J) + OMEGA/4.0*( PSI(I-1,J) + PSI(I+1,J) + PSI(I,J-1) + PSI(I,J+1) );
            ERROR = ERROR + abs(PSI(I,J)-PSIPRV);
        end
    end
end


% ***** PLOT REPRESENTATIVE STREAMLINES *****

% [Q,P] = contour(PSI');
% clabel(Q);
% title([['Pattern of flow through a chamber: OMEGA = ' num2str(OMEGA)] [' and number of iterations = ' num2str(ITER)]])
% xlabel('X - AXIS')
% ylabel('Y - AXIS')

X(1) = 0.0;
for I=2:M
    X(I) = X(I-1) + H;
end

Y(1) = 0.0;
for J=2:N
    Y(J) = Y(J-1) + H;
end

hold on;
for K=1:length(CONST)
    [XX,YY,KMAX] = SEARCH3(X,Y,PSI,CONST(K));
    plot(XX(1:KMAX),YY(1:KMAX),SYMBOL(K))
end
hold off;

title([['Pattern of flow through a chamber: OMEGA = ' num2str(OMEGA)] [' and number of iterations = ' num2str(ITER)]])
xlabel('X - AXIS')
ylabel('Y - AXIS')
legend(['PSI=' num2str(CONST(1))],['PSI=' num2str(CONST(2))],['PSI=' num2str(CONST(3))],['PSI=' num2str(CONST(4))],['PSI=' num2str(CONST(5))],['PSI=' num2str(CONST(6))],['PSI=' num2str(CONST(7))],['PSI=' num2str(CONST(8))],['PSI=' num2str(CONST(9))],['PSI=' num2str(CONST(10))],['PSI=' num2str(CONST(11))])
