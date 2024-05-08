function [XX, YY, KMAX] = SEARCH(X, Y, PSI, PSIA);

% ***** SEARCHING THROUGH VERTICAL GRID LINES FOR POINTS (XX,YY) ON A 
% STREAMLINE ASSUMING THE VALUE PSIA. PSI CONTAINS VALUES OF STREAM
% FUNCTION AT GRID POINTS (X,Y). KMAX IS THE TOTAL NUMBER OF POINTS THUS 
% FOUND *****

global M N W DY;

% ***** EXAMINE EVERY INTERVAL ON EACH OF THE M VERTICAL GRID LINES. A
% POINT (XX,YY) IS LOCATED AT A GRID POINT IF THE LOCAL PSI DIFFERS FROM
% THE DESIRED VALUE BY 0.00001. OTHERWISE WE LOOK FOR THE INTERVALS ACROSS
% WHICH P*Q<=0, AND THEN COMPUTE THE COORDINATES (XX,YY) WITHIN THOSE
% INTERVALS *****

XX = zeros(W,1);
YY = zeros(W,1);

K = 0;
for I=1:M

    J = 1;
    while (J <= N)

        P = PSI(I,J) - PSIA;
        if ( abs(P) <= 0.0001 )
            K = K + 1;
            XX(K) = X(I);
            YY(K) = Y(J);
        else
            while (J < N)
    
                J = J + 1;
                Q = PSI(I,J) - PSIA;
                if ( P*Q > 0.0 )
                    P = Q;
                else
                    K = K + 1;
                    XX(K) = X(I);
                    YY(K) = Y(J) - DY*abs(Q)/(abs(P)+abs(Q));
                    P = PSI(I,J) - PSIA;
                end

            end
        end

        J = J + 1;
    end

end

KMAX = K;
