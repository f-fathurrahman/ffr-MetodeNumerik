%                       ***** PROGRAM 2.1 *****                            
% SOLVING A BOUNDARY-VALUE PROBLEM CONCERNING THE AXISYMMETRIC FLOW CAUSED
% BY A SOURCE DISTRIBUTION AND A LINE SINK AT ORINGIN
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clear all; close all; clc
format short e;

RMIN = 1.0;
RMAX = 4.0;
INTRVL = 30;


% ***** USE THE ARITHMETIC FUNCTIONS TO REPRESENT (2.3.10) *****

A = @(X)1./X;
B = @(X)0.*X;
D = @(X)9.*X;


% ***** N AND H DEFINED IN FIG. 2.2.1 ARE COMPUTED *****

N = INTRVL - 1;
H = (RMAX-RMIN) / INTRVL;


% ***** COMPUTE R(I), C(I,J), AND EXACT VALUES OF PHI AND U *****

RADIUS = RMIN + H;

for I=1:N
    R(I) = RADIUS;
    C(I,1) = 1. - H/2.*feval(A,R(I));
    C(I,2) = -2. + H*H*feval(B,R(I));
    C(I,3) = 1. + H/2.*feval(A,R(I));
    C(I,4) = H*H*feval(D,R(I));
    EXPHI(I) = R(I)^3 - 24.*log(R(I));
    EXU(I) = 3.*R(I)^2 - 24./R(I);
    RADIUS = RADIUS + H;   
end


% ***** SPECIFY BOUNDARY CONDITIONS AND MODIFY SOME OF THE COEFFICIENTS 
% ACCORDING TO (2.2.12) *****

PHI0 = 1.0;
PHINP1 = 64. - 24.*log(4.);
C(1,4) = C(1,4) - C(1,1)*PHI0;
C(1,1) = 0.;
C(N,4) = C(N,4) - C(N,3)*PHINP1;
C(N,3) = 0.;


% ***** COMPUTE PHI AND U AT INTERIOR POINTS *****

[C,PHI] = TRID(C,N);
U(1) = (PHI(2)-PHI0) / (2.*H);
NM1 = N-1;
for J=2:NM1
    U(J) = (PHI(J+1)-PHI(J-1)) / (2.*H);
end
U(N) = (PHINP1-PHI(NM1)) / (2.*H);


% ***** PRINT THE RESULT IN COMPARISON WITH EXACT SOLUTION *****

disp('       I           R(I)        PHI(I)     EXPHI(I)        U(I)        EXU(I)');
for I=1:N
    disp([I R(I) PHI(I) EXPHI(I) U(I) EXU(I)]);
end

plot(R,U)
hold on
plot(R,PHI,'--')
xlabel('R '),ylabel ('U, \phi')
legend(' U, ---\phi')
title('Radial jet flow')
