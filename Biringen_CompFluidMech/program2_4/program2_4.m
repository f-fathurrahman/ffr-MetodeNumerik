%                       ***** PROGRAM 2.4 *****                            

% TRANSFORMATION FROM THE FLOW PAST A CIRCULAR CYLINDER TO THAT AROUND A 
% JOUKOWSKI AIRFOIL THROUGH CONFORMAL MAPPING
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clf;
clear all;
format short e;

global M N W DY A B C ZP1;

XMAX = 3.;
XMIN = -3.;
YMAX = 3.;
YMIN = -3.;
XMAX1 = 3.;
XMIN1 = -3.;
YMAX1 = 3.;
YMIN1 = -3.;
DX1 = 0.05;
DY1 = 0.05;
M = (XMAX1-XMIN1)/DX1 + 1.;
N = (YMAX1-YMIN1)/DY1 + 1.;
W = 200;

X = zeros(W,1);
Y = zeros(W,1);
X1 = zeros(M,1);
Y1 = zeros(N,1);
PSI = zeros(M,N);
CONST = zeros(10,1);
PHI = zeros(W,1);
CP = zeros(W,1);


% ***** COMPUTE XP1 and ZP1 FROM GIVEN VALUES OF A, B, U AND XP1 *****

A = 1.0;
B = 0.8;
C = 0.99;
U = 1.0;
CI = i;
YP1 = 0.199;
XP1 = B - sqrt(A^2-YP1^2);
ZP1 = XP1 + YP1*i;


% ***** DEFINE GRID POINTS IN Z1-PLANE AND EVALUATE STREAM FUNCTION THERE *****

DY = DY1;
X1(1) = XMIN1;
for I=2:M
    X1(I) = X1(I-1) + DX1;
end
Y1(1) = YMIN1;
for J=2:N
    Y1(J) = Y1(J-1) + DY1;
end

for I=1:M
    for J=1:N
        Z1 = X1(I) + Y1(J)*i;
        Z2 = Z1 - ZP1;
        PSI(I,J) = U*imag( Z2 + A^2/Z2 + CI*2.*YP1*log(Z2/A));
     end
end


% ***** SEARCH FOR POINTS ON THE FIRST STREAMLINE IN Z1-PLANE *****

% [Q,P] = contour(PSI');
% clabel(Q);
% title('Stream function plot.')
% xlabel('X - AXIS')
% ylabel('Y - AXIS')

CONST(1) = -2.;
[X,Y,KMAX] = SEARCH(X1,Y1,PSI,CONST(1));
[XX, YY, NPOINT] = MAPPNG(X, Y, KMAX, XMIN, XMAX, YMIN, YMAX);
plot(XX(1:NPOINT),YY(1:NPOINT),'.','markersize',8)
axis equal
title('Flow around a Joukowski airfoil')
xlabel('X - AXIS')
ylabel('Y - AXIS')

hold on;
for K=2:length(CONST)
    CONST(K) = CONST(K-1) + 0.5;
    [X,Y,KMAX] = SEARCH(X1,Y1,PSI,CONST(K));
    [XX, YY, NPOINT] = MAPPNG(X, Y, KMAX, XMIN, XMAX, YMIN, YMAX);
    plot(XX(1:NPOINT),YY(1:NPOINT),'.')
end
hold off;


% ***** SPECIFY THE ANGLE PHI FOR 100 POINTS ON THE CIRCLE WHI *****

DPHI = pi/100;
for L=1:W
    if (L > 1)
        PHI(L) = PHI(L-1) + DPHI;
    else
        PHI(L) = 0.0;
    end

    Z2 = A * exp(CI*PHI(L));
    Z1 = Z2 + ZP1;
    Z = Z1 + B^2/Z1;
    X(L) = real(Z);
    Y(L) = imag(Z);
    V = U*abs( (1.-(A/Z2)^2+CI*2.*YP1/Z2)/(1.-(B/Z1)^2) );
    CP(L) = 1. - (V/U)^2;
end

figure;
plot(X,CP,'.','markersize',8)
axis equal
title('Pressure distribution around the airfoil')
xlabel('X - AXIS')
ylabel('CP : PRESSURE COEFFICIENT')



