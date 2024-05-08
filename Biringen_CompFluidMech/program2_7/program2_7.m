%                       ***** PROGRAM 2.7 *****                            

% CHANNEL FLOW PAST A CIRCULAR CYLINDER, AN EXAMPLE ON CURVED SURFACE AND
% THE DERIVATIVE BOUNDARY CONDITION. PLOTS OBTAINED USING THE SEARCH PROGRAM
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc; clear all; close all;
clear all;
format short e;

global M N H JFIRST;

M = 25;
N = 9;

SYMBOL = zeros(11,1);
CONST = zeros(11,1);
JFIRST = zeros(M,1);
ID = zeros(M,N);
PSI = zeros(M,N);
PSIF = zeros(M,N);
A = zeros(M,N);
B = zeros(M,N);
C = zeros(M,N);
D = zeros(M,N);
X = zeros(M,1);
Y = zeros(N,1);


% ***** ASSIGN INPUT DATA *****

H = 0.25;
SYMBOL = '+osdv^x*ph.';
CONST = [0.0; 0.1; 0.2; 0.3; 0.4; 0.5; 0.6; 0.7; 0.8; 0.9; 1.0];


% ***** DEFINE THE VALUE JAT THE LOWEST INTERIOR GRID POINT ALONG EACH 
% VERTICAL GRID LINE *****

for I=1:M
    JFIRST(I) = 2;
    if (I >= 11 && I <= 15)
        JFIRST(I) = 5;
    end
end

JFIRST(10) = 4;
JFIRST(13) = 6;
JFIRST(16) = 4;


% ***** ASSIGN PSI AT GRID POINTS ON OR INSIDE SOLID BOUNDARIES *****

for I=1:M
    J1 = JFIRST(I) - 1;
    for J=1:J1
        PSI(I,J) = 0.0;
    end
end

for I=1:M
    PSI(I,N) = 1.0;
end

for J=2:N-1
    PSI(1,J) = 1.0;
end


% ***** ASSIGN ID NUMBERS TO ALL INTERIOR POINTS *****

for I=2:M
    JA = JFIRST(I);
    for J=JA:N-1
        ID(I,J) = 0;
        if (I == M)
            ID(I,J) = -1;
        end
    end
end

for I=9:17
    J = JFIRST(I);
    ID(I,J) = 1;
end

ID(9,3) = 1;
ID(13,6) = 0;
ID(17,3) = 1;


% ***** COMPUTE A, B, C AND D FOR THE POINTS AT WHICH (2.11.8) WILL BE 
% USED *****

for I=9:17
    for J=2:5
        if (ID(I,J) == 1)
            A(I,J) = H;
            B(I,J) = H;
            C(I,J) = H;
            D(I,J) = H;
            if (I == 9 || I == 10)
                B(I,J) = (13-I)*H-sqrt(1.-((J-1)*H)^2);
            end
            if (I == 16 || I == 17)
                A(I,J) = (I-13)*H-sqrt(1.-((J-1)*H)^2);
            end
            if (I >= 10 && I <= 16)
                C(I,J) = (J-1)*H-sqrt(1.-((13-I)*H)^2);
            end
        end
    end
end


% ***** ASSIGN GUESSED VALUES TO PSI *****

for I=2:M
    JA = JFIRST(I);
    for J=JA:N-1
        PSI(I,J) = 0.5;
    end
end


% ***** ITERATE UNTIL SUM OF ERROS IS LESS THAN OR EQUAL TO 0.0005 *****

ITER = 0;
ERROR = 1.0;
while (ERROR > 0.0005)
    ITER = ITER + 1;
    ERROR = 0.0;
    for I=2:M
        JA = JFIRST(I);
        for J=JA:N-1
            PSIPRV = PSI(I,J);
            if (ID(I,J) == -1)
%(2.11.13) IS USED TO INCORPORATE DERIVATIVE BOUNDARY CONDITION AT THE EXIT:
                PSI(I,J) = 0.25 * ( 2.*PSI(I-1,J) + PSI(I,J-1) + PSI(I,J+1) );                
            elseif (ID(I,J) == 0)
%AT GRID POINTS AWAY FROM THE CURVED BOUNDARY (2.11.6) IS USED:
                PSI(I,J) = 0.25 * ( PSI(I-1,J) +  PSI(I+1,J) + PSI(I,J-1) + PSI(I,J+1) );
            elseif (ID(I,J) == 1)
%AT POINTS NEIGHBORING CURVED BOUNDARY (2.11.8) IS USED:
                PSI(I,J) = ( PSI(I-1,J)/(A(I,J)*(A(I,J)+B(I,J))) +  PSI(I+1,J)/(B(I,J)*(A(I,J)+B(I,J))) + PSI(I,J-1)/(C(I,J)*(C(I,J)+D(I,J))) + PSI(I,J+1)/(D(I,J)*(C(I,J)+D(I,J))) ) / ( 1./(A(I,J)*B(I,J)) + 1./(C(I,J)*D(I,J)) );
            end
            ERROR = ERROR + abs(PSI(I,J)-PSIPRV);
        end
    end
end


% ***** PLOT REPRESENTATIVE STREAMLINES *****

% [Q,P] = contour(PSI');
% clabel(Q);
% title(['Chanel flow past a circular cylinder: number of iterations = ' num2str(ITER)])
% xlabel('X - AXIS')
% ylabel('Y - AXIS')

% figure;
% for I=1:M
%    for J=N:-1:1
%       PSIF(I,N-J+1)=PSI(I,J); 
%    end
% end
% [QF,PF] = contour3(PSIF');
% clabel(QF);
% title('Stream function plot.')
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
    [XX,YY,KMAX] = SEARCH4(X,Y,PSI,CONST(K));
    plot(XX(1:KMAX),YY(1:KMAX),SYMBOL(K))
end
hold off;

title(['Chanel flow past a circular cylinder: number of iterations = ' num2str(ITER)])
xlabel('X - AXIS')
ylabel('Y - AXIS')
legend(['PSI=' num2str(CONST(1))],['PSI=' num2str(CONST(2))],['PSI=' num2str(CONST(3))],['PSI=' num2str(CONST(4))],['PSI=' num2str(CONST(5))],['PSI=' num2str(CONST(6))],['PSI=' num2str(CONST(7))],['PSI=' num2str(CONST(8))],['PSI=' num2str(CONST(9))],['PSI=' num2str(CONST(10))],['PSI=' num2str(CONST(11))])

hold on;
for K=1:length(CONST)
    [XX,YY,KMAX] = SEARCH4(X,Y,PSI,CONST(K));
    plot(XX(1:KMAX),-YY(1:KMAX),SYMBOL(K))
end
hold off;
