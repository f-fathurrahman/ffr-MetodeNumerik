%                       ***** PROGRAM 2.3 *****                            

% APPROXIMATE THE FLOW OF A UNIFORM STREAM PAST A SPRERE BY THE USE OF  
% VON KARMAN'S METHOD 
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clear all; close all; clc
format short e;

CONST = zeros(10,1);
CP = zeros(9,1);
CPEX = zeros(9,1);
D = zeros(9,11);
DD = zeros(10,10);
PHI = zeros(9,1);
Q = zeros(10,1);
R = zeros(9,1);
S = zeros(10,1);
V = zeros(9,1);
VEX = zeros(9,1);
Z = zeros(9,1);
ZZ = zeros(11,1);
A = 1.0;
U = 1.0;


% ***** COMPUTE THE POSITIONS ON THE Z-AXIS OF THE END POINTS OF SOURCE 
% SEGMENTS, WHICH ARE OF THE SAME LENGTH OF 0.16A *****

ZZ(1) = -0.8 * A;
for I=2:11
    ZZ(I) = ZZ(I-1) + 0.16*A;
    S(I-1) = 0.16*A;
end


% ***** COMPUTE THE COORDINATES OF THE SURFACE POINTS ACCORDING TO (2.5.21
% AND 22
% *****

for I=1:9
    PHI(I) = ( 10. - I ) * pi / 10.;
    R(I) = A * sin(PHI(I));
    Z(I) = A * cos(PHI(I));
end


% ***** COMPUTE D(I,J) AND DD(I,J) ACCORDING TO (2.6.8) *****

for I=1:9
    for J=1:11
        D(I,J) = sqrt( R(I)^2+(Z(I)-ZZ(J))^2 );    
    end
end

for I=1:9
    for J=1:10
        DD(I,J) = D(I,J) - D(I,J+1);
    end
end


% ***** IN (2.6.11), CALL THE COEFFICIENT MATRIX DD(I,J) AND CALL THE
% COLUMN MATRIX ON THE RIGHT CONST(I). PRINT OUT THESE TWO MATRICES AFTER
% ALL ELEMENTS HAVE BEEN COMPUTED *****

disp('THE SQUARE MATRIX ON THE LEFT AND THE COLUMN MATRIX ON THE RIGHT OF (2.6.11) ARE :');
for J=1:10
    DD(10,J) = S(J);
end
for I=1:9
    CONST(I) = R(I)^2 / 2.0; 
end
CONST(10) = 0.0;
for I=1:10
    disp([ [DD(I,1:10)] CONST(I) ]);
end

disp('SOURCE STRENGTHS Q(I) ARE, FOR I = 1,2,...,10 :');
%Q = DD \ CONST;
Q = inv(DD) * CONST;
disp(Q');

% ***** COMPUTE VELOCITY AND PRESSURE COEFFICIENT AT EACH SURFACE POINT
% USING (2.6.12, 13 AND 14) AND THE EXACT VALUES USING (2.6.20 AND 21).
% PRINT OUT THE RESULTS FOR COMPARISON *****

disp('COMPARISON OF NUMERICAL RESULT WITH THE EXACT AT 9 SURFACE POINTS');
disp('      R(I)        Z(I)         V(I)         VEX(I)       CP(I)       CPEX(I)')
for I=1:9
    UR = 0.0;
    UZ = 1.0;
    for J=1:10
        UR = UR + Q(J)/R(I)*((Z(I)-ZZ(J+1))/D(I,J+1)-(Z(I)-ZZ(J))/D(I,J));
        UZ = UZ + Q(J)*(1./D(I,J+1)-1./D(I,J));
    end
    V(I) = U * sqrt(UR^2+UZ^2);
    VEX(I) = 1.5*U*R(I)/A;
    CP(I) = 1. - (V(I)/U)^2;
    CPEX(I) = 1. - (VEX(I)/U)^2;
    disp([R(I) Z(I) V(I) VEX(I) CP(I) CPEX(I)]);
end

