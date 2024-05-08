%                       ***** PROGRAM 2.5 *****                            

% SOLVING POISSON EQUATION FOR THE FLOW AROUND A VORTEX BOUNDED WITHIN A 
% RECTANGULAR REGION USING LIEBMANN'S ITERATIVE FORMULA 
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clear all; close all; clc; clf;
clear all;
format short e;

global M N H;

XMAX = 3.;
XMIN = -3.;
YMAX = 2.;
YMIN = -2.;
H = 0.25;
M = (XMAX-XMIN)/H + 1.;
N = (YMAX-YMIN)/H + 1.;

X = zeros(M,1);
Y = zeros(N,1);
XX = zeros(100,1);
YY = zeros(100,1);
CONST = zeros(5,1);
PLTC = zeros(5,1);
PSI = zeros(M,N);
PSIPRV = zeros(M,N);
ZETA = zeros(M,N);


% ***** ASSIGN INPUT DATA *****

ZETA0 = 100.;
X0 = 1.;
Y0 = 1.;
ITMAX = 200;
ERRMAX = 0.001;
CONST = [0.05; 0.2; 0.5; 1.0; 1.5]
PLTC = 'osdv^+';


% ***** COMPUTE COORDINATES FOR GRID POINTS *****

X(1) = XMIN;
for I=2:M
    X(I) = X(I-1) + H;
end
Y(1) = YMIN;
for J=2:N
    Y(J) = Y(J-1) + H;
end


% ***** ASSIGN BOUNDARY VALUES TO PSI *****

for I=1:M
    PSI(I,1) = 0.;
    PSI(I,N) = 0.;
end
    
for J=2:N-1
    PSI(1,J) = 0.;
    PSI(M,J) = 0.;
end


% ***** CONVERT COORDINATES (X0,Y0) OF THE VORTEX TO (I0,J0) IN INDEX 
% NOTATION. THEN ASSIGN VALUES TO ZETA WHICH HAS A NONVANISHING VALUE OF
% ZETA0 AT (X0,Y0) *****

I0 = (X0-XMIN)/H + 1.;
J0 = (Y0-YMIN)/H + 1.;
for I=1:M
    for J=1:N
        ZETA(I,J) = 0.;
    end
end
ZETA(I0,J0) = ZETA0;


% ***** ASSIGN GUESSED VALUES TO PSI *****

for I=2:M-1
    for J=2:N-1
        PSI(I,J) = 0.5;
    end
end


% ***** LET ITERATION COUNTER START FROM ZERO *****

ITER = 0;
ERROR = 0.1;


% ***** START AN ITERATION BY INCREASING ITERATION COUNTER BY ONE AND 
% SETTING ERROR INITIALLY TO ZERO, BEFORE APPLYING LIEBMANN'S FORMULA, THE 
% LOCAL VALUES OF PSI ARE STORED IN THE ARRAY PSIPRV. ABSOLUTE DIFFERENCES 
% BETWEEN TWO CONSECUTIVE APPROXIMATIONS TO PSI AT INDIVIDUAL INTERIOR 
% POINTS ARE SUMMED AND STORED IN ERROR *****

for I=1:M
    PSIPRV(I,1) = PSI(I,1);
    PSIPRV(I,N) = PSI(I,N);
end

for J=2:N-1
    PSIPRV(1,J) = PSI(1,J);
    PSIPRV(M,J) = PSI(M,J);
end


% ***** PRINT PSI IF FINAL VALUE OF ERROR IS LESS THAN OF EQUAL TO ERRMAX, 
% OR IF THE NUMBER OF ITERATIONS HAS REACHED THE VALUE ITMAX, OTHERWISE GO 
% BACK FOR ANOTHER ITERATION *****

while (ERROR > ERRMAX && ITER < ITMAX)
    ITER = ITER + 1;
    ERROR = 0.;
    for I=2:M-1
        for J=2:N-1
            PSIPRV(I,J) = PSI(I,J);
            PSI(I,J) = (PSI(I-1,J) + PSI(I+1,J) + PSI(I,J-1) + PSI(I,J+1) + H^2*ZETA(I,J))/4.0;
            ERROR = ERROR + abs(PSI(I,J)-PSIPRV(I,J));
        end
    end
end

disp([[[['M = ' num2str(M)] ['     N = ' num2str(N)]] ['    NUMBER OF ITERATIONS = ' num2str(ITER)]] ['     SUM OF ERRORS = ' num2str(ERROR)]])

J = 0;
for K=1:2:N
    L = N - K + 1;
    J = J + 1;
    I = 0;
    for P=1:3:M
        I = I + 1;
        PSIMOD(I,J) = PSI(P,L);
    end
end

disp('THE RESULT FOR PSI(I,J):');
disp(PSIMOD');

% ***** PRINT PSIPRV TO MAKE SURE THAT THE RESULT IS SATISFACTORY AT EVERY 
% GRID POINT, THEN SKIP TO THE NEXT PAGE *****

J = 0;
for K=1:2:N
    L = N - K + 1;
    J = J + 1;
    I = 0;
    for P=1:3:M
        I = I + 1;
        PSIPRVMOD(I,J) = PSIPRV(P,L);
    end
end

disp('THE RESULT FOR PSI(I,J):');
disp(PSIPRVMOD');


% ***** PLOT FIVE REPRESENTATIVE STREAMLINES *****

% [Q,P] = contour(PSI');
% clabel(Q);
% title('Streamlines of a vortex bounded by a rectangular wall.')
% xlabel('X - AXIS')
% ylabel('Y - AXIS')

hold on;
for I=1:length(CONST)
    [XX,YY,KMAX] = SEARCH2(X,Y,PSI,CONST(I));
    plot(XX(1:KMAX),YY(1:KMAX),PLTC(I))
end
plot(X0,Y0,PLTC(I+1))
hold off;

%title('Streamlines of a vortex bounded by a rectangular wall.')
xlabel('X - AXIS')
ylabel('Y - AXIS')
legend(['PSI=' num2str(CONST(1))],['PSI=' num2str(CONST(2))],['PSI=' num2str(CONST(3))],['PSI=' num2str(CONST(4))],['PSI=' num2str(CONST(5))])
