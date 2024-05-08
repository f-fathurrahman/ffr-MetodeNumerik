%                       ***** PROGRAM 2.9 *****                            

% PROPAGATION OF A FINITE-AMPLITUDE WAVE IN A TUBE HAVING A CLOSED LEFT END
% AND AN OPEN RIGHT END
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clf;
clear all;
format short e;

M = 101;
JMAX = 55;

A = zeros(M,JMAX+1);
U = zeros(M,JMAX+1);


% ***** ASSIGN INPUT DATA *****

A0 = 340.;
H = 0.02;
GAMMA = 1.4;
NT = 5;
TAU = 0.5 * H/A0;
RATIO = TAU / H;


% ***** ASSIGN INITIAL CONDITION TO U. FOR A RIGHTWARD-TRAVELING WAVE, 
% CHOOSE THE UPPER SIGN IN (2.12.15) TO COMPUTE A *****

AMPLTD = A0 / 2.0;
COEFF = (GAMMA-1.0) / 2.0;
for I=1:M
    U(I,1) = 0.0;
    if (I >= 1 && I <= 13)
        U(I,1) = AMPLTD * (I-1)/12;
    end
    if (I > 13 && I <= 39)
        U(I,1) = AMPLTD * (26-I)/13;
    end
    if (I > 39 && I <= 51)
        U(I,1) = AMPLTD * (I-51)/12;
    end
    A(I,1) = A0 + COEFF*U(I,1);
end


% ***** ASSIGN BOUNDARY CONDITION FOR U ON THE CLOSED LEFT END AND THAT FOR
% A ON THE OPEN RIGHT END*****

for J=1:JMAX
    U(1,J+1) = 0.0;
    A(M,J+1) = A0;
end


% ***** COMPUTE SOLUTION USING THE METHOD OF COURANT, ISAACSON AND REES. 
% SPECIAL FORMULAE ARE USED AT THE BOUNDARIES WHERE I=1 AND I=M, STABILITY
% CONDITION FOR THE NUMERICAL SCHEME IS CHECKED AT EVERY GRID POINT *****

JLAST = 0;
for J=1:JMAX
    UPA = RATIO * (U(1,J)+A(1,J));
    UMA = RATIO * (U(1,J)-A(1,J));
    if (abs(UPA) <= 1. || abs(UMA) <= 1.)
        UB = U(1,J) - UMA*(U(2,J)-U(1,J));
        AB = A(1,J) - UMA*(A(2,J)-A(1,J));
        A(1,J+1) = AB - COEFF*UB;
        for I=2:M-1
            UPA = RATIO * (U(I,J)+A(I,J));
            UMA = RATIO * (U(I,J)-A(I,J));
            if (abs(UPA) <= 1. || abs(UMA) <= 1.)
                UA = U(I,J) + UPA*(U(I-1,J)-U(I,J));
                AA = A(I,J) + UPA*(A(I-1,J)-A(I,J));
                UB = U(I,J) - UMA*(U(I+1,J)-U(I,J));
                AB = A(I,J) - UMA*(A(I+1,J)-A(I,J));
                U(I,J+1) = 0.5 * ( (UA+UB) + (AA-AB)/COEFF );
                A(I,J+1) = 0.5 * ( COEFF*(UA-UB) + (AA+AB) );
            else
                JLAST = J;
                break;
            end
        end
        if (JLAST == 0)
            UPA = RATIO * (U(M,J)+A(M,J));
            UMA = RATIO * (U(M,J)-A(M,J));
            if (abs(UPA) <= 1. || abs(UMA) <= 1.)
                UA = U(M,J) + UPA*(U(M-1,J)-U(M,J));
                AA = A(M,J) + UPA*(A(M-1,J)-A(M,J));
                U(M,J+1) = UA + (AA-A0)/COEFF;
            else
                JLAST = J;
                break;
            end
        end
    else
        JLAST = J;
        break;
    end
    if (JLAST ~= 0)
        break;
    end
end


% ***** COMPUTE COORDINATES FOR GRID POINTS AND PLOT SOLUTION *****

X(1) = 0.0;
for I=2:M
    X(I) = X(I-1) + H;
end

if (JLAST == 0)
    JLAST = JMAX+1;
end
NN=1;
figure;
for J=1:5:26
    subplot(6,1,NN)
    plot(X,U(:,J))
    NN=NN+1;
    set(gca,'YTick',[-150 0 150])
    set(gca,'XTick',[0 2])
    axis([0,2,-150,150])
    grid on
    title([['Time = ' int2str(J-1)] '*TAU = ' num2str((J-1)*TAU) ' seconds'])
end

NN=1;
figure;
for J=31:5:56
    subplot(6,1,NN)
    plot(X,U(:,J))
    NN=NN+1;
    set(gca,'YTick',[-150 0 150])
    set(gca,'XTick',[0 2])
    axis([0,2,-150,150])
    grid on
    title([['Time = ' int2str(J-1)] '*tau = ' num2str((J-1)*TAU) ' seconds'])
end
