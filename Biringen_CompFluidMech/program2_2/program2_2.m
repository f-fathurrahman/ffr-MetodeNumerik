%                       ***** PROGRAM 2.2 *****                            

% PLOTTING THE FLOW RESULTING FROM SUPERPOSITION OF A UNIFORM STREAM AND 
% FIVE DOUBLETS 
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clf;
clear all;
format short e;

global M N W DY;

% ***** ASSIGN INPUT DATA *****

XMAX = 3.;
XMIN = -3.;
YMAX = 2.;
YMIN = -2.;
DX = 0.2;
DY = 0.2;
X = XMIN:DX:XMAX;
Y = YMIN:DY:YMAX;
M = length(X);
N = length(Y);
W = 100;

XX=zeros(W,1);
YY=zeros(W,1);
PSI=zeros(M,N);

CONST = [0., -0.5, -1.0, 0.5, 1.0];
PLTC = ['ks';'ko';'kd';'kv';'k^'];


% ***** EVALUATE THE VALUES OF STREAM FUNCTION AT THE GRID POINTS *****

for I=1:M
    for J=1:N
        [PSI(I,J)] = F(X(I),Y(J));
    end
end


% ***** SEARCH FOR POINTS ON 5 STREAMLINES ALONG WHICH THE STREAM FUNCTION 
% HAS THE VALUES STORED IN CONST. THEIR COORDINATES ARE PRINTED AND ARE 
% REPRESENTED ON THE GRAPH BY 5 DIFFERENT PLOTTING CHARACTERS STORED IN 
% PLTC *****

% [Q,P] = contour(PSI');
% clabel(Q);
% title('Stream function plot.')
% xlabel('X - AXIS')
% ylabel('Y - AXIS')

hold on;
for L=1:length(CONST)
    PSIA = CONST(L);
    [XX,YY,KMAX] = SEARCH2(X,Y,PSI,PSIA);
    disp(['THE POINTS (X,Y) ON THE STREAMLINE PSI=' num2str(PSIA)]);
    disp([XX(1:KMAX)' YY(1:KMAX)'])
    plot(XX(1:KMAX),YY(1:KMAX),PLTC(L,:),'MarkerSize',4)    
end
hold off;

legend(['PSI=' num2str(CONST(1))],['PSI=' num2str(CONST(2))],['PSI=' num2str(CONST(3))],['PSI=' num2str(CONST(4))],['PSI=' num2str(CONST(5))]);
