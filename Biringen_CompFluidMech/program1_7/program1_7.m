%                       ***** PROGRAM 1.7 *****                            

% STUDYING THE ROLL-UP OF A VORTEX SHEET
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clf;
clear all;
format short e;

global I M G X Y;

M = 40;
NMAX = 500;
NFREQ = 5;

X = zeros(M,1);
Y = zeros(M,1);
XNEW = zeros(M/2,1);
YNEW = zeros(M/2,1);
U = zeros(M/2,1);
V = zeros(M/2,1);
UNEW = zeros(M/2,1);
VNEW = zeros(M/2,1);
CONST = zeros(7,1);

CONST = [0.005; 0.025; 0.050; 0.075; 0.10; 0.15; 0.20; 0.35];

% ***** SPECIFY DIMENSIONLESS STRENGTH AND INITIAL COORDINATES OF VORTICES
% *****

NDT = 0;
T = 0.0;
DX = 1. / M;
X(1) = -1. + DX;
for I=2:M
    X(I) = X(I-1) + 2.*DX;
end
for I=1:M/2
    Y(I) = 0.0;
    G(I) = (1.-(X(I)+DX)^2)^0.5 - (1.-(X(I)-DX)^2)^0.5;
end
for I=M/2+1:M
    Y(I) = 0.0;
    G(I) = -G(M+1-I);
end

% ***** COMPUTE THE POSITION OF EACH VORTEX ON THE LEFT HALF OF THE SHEET
% AT TIME DT LATER *****

DT = 1. / (25.*M);
K = 1;
while (NDT <= NMAX)

    NDT = NDT + 1;
    T = T + DT;
    for I=1:M/2
        [TSOL,SOL] = ode45(@UV, [T-DT, T], [X(I), U(I), Y(I), V(I)]);
        XNEW(I)=SOL(end,1);
        UNEW(I)=SOL(end,2);
        YNEW(I)=SOL(end,3);
        VNEW(I)=SOL(end,4);
    end

    
% ***** ASSIGN VALUES OF XNEW AND YNEW TO X AND Y, RESPECTIVELY *****

    for I=1:M/2
        X(I) = XNEW(I);
        Y(I) = YNEW(I);
        U(I) = UNEW(I);
        V(I) = VNEW(I);
    end


% ***** VORTICES ON THE RIGHT HALF OF THE SHEET ARE THE MIRROR IMAGES OF 
% THOSE ON THE LEFT *****    

    for I=M/2+1:M
        X(I) = -X(M+1-I);
        Y(I) = Y(M+1-I);
    end
    
    if (mod(NDT,NFREQ) == 0)
        disp(['NDT = ' num2str(NDT)]);
        disp('      X            Y          ');
        disp('    ------       ------       ');
        disp([X(M/2+1:M) Y(M/2+1:M)]);
    end

    if (K <= length(CONST))
        if (T >= CONST(K))
            plot(X(M/2+1:M),Y(M/2+1:M),'.')
            title([['Evolution of an initially flat vortex sheet T = ' num2str(T)] '.'])
            xlabel('X - AXIS')
            ylabel('Y - AXIS')
            K = K + 1;
            if (K <= length(CONST))
                figure;
            end
        end
    end

end    
    

