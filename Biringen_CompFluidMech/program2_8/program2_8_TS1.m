%                       ***** PROGRAM 2.8 *****                            

% LEAPFROG METHOD IS USED TO STUDY THE PROPAGATION OF A SINUSOIDAL WAVELET
% IN A TUBE WITH LEFT END CLOSED AND RIGHT END OPEN TO THE ATMOSPHERE
% Reference: Computational Fluid Mechanics by Example (Biringen and Chow; Wiley 2011)

clc;
clf;
clear all;
format short e;

M = 51;
JMAX = 116;


F = zeros(M,1);
G = zeros(M,1);
U = zeros(M,JMAX);



% ***** ASSIGN INPUT DATA *****

A = 340.;
H = 0.02;
NT = 5;
TAU = H / A;


% ***** DEFINE FUNCTIONS F AND G *****

for I=1:M
    F(I) = 0.0;
    G(I) = 0.0;
    if (I >= 11 && I <= 31)
        F(I) = sin( 5.*pi*((I-1)*H-0.2) );
    end
end


% ***** ASSIGN INITIAL CONDITION TO U *****

for I=1:M
    U(I,1) = F(I);
end


% ***** ASSIGN BOUNDARY CONDITION AT THE CLOSED LEFT END *****

for J=2:JMAX
    U(1,J) = 0.0;
end


% ***** COMPUTE SOLUTION AT J=2 USING STARTING FORMULAE (2.12.18) AND 
% (2.13.9) *****

for I=2:M-1
    U(I,2) = (F(I-1)+F(I+1))/2.0 + TAU*G(I);
end
U(M,2) = F(M-1) + TAU*G(M);


% ***** THEN COMPUTE SOLUTION AT SUBSEQUENT TIME STEPS USING (2.12.13) AND
% (2.12.13). DERIVATIVE BOUNDARY CONDITION IS INCORPORATED IN THE LATTER 
% FORMULA *****

for J=2:JMAX-1
    for I=2:M-1
        U(I,J+1) = U(I-1,J) + U(I+1,J) - U(I,J-1);
    end
    U(M,J+1) = 2.0*U(M-1,J) - U(M,J-1);
end

for J=1:JMAX
    for I=1:M
        if (abs(U(I,J))<1e-10)
            U(I,J) = 0.0;
        end
    end
end


% ***** COMPUTE COORDINATES FOR GRID POINTS AND PLOT SOLUTION *****

X(1) = 0.0;
for I=2:M
    X(I) = X(I-1) + H;
end
%zsz=zeros(size(X));
NN=1
figure
for J=1:5:16
        subplot(4,1,NN)
        plot(X,U(:,J));
        %plot(X,zsz,'k')
        NN=NN+1
     
    axis([0,1,-1,1])
    
    
    set(gca,'YTick',[-1 0 1])
    set(gca,'XTick',[0  1])
    grid on
    title([['Time = ' int2str(J-1)] '*tau = ' num2str((J-1)*TAU) ' seconds'])
end
figure;
    NNN=1
    for J=21:5:36
        
        subplot(4,1,NNN)
        plot(X,U(:,J))
        NNN=NNN+1
        set(gca,'YTick',[-1 0 1])
    set(gca,'XTick',[0 0.5 1])
    grid on;
        axis([0,1,-1,1])
    %xlabel('X')
    %ylabel('Magnitude')
    %title([['Time = ' num2str((J-1)*TAU)] ' seconds'])
    title([['Time = ' int2str(J-1)] '*tau = ' num2str((J-1)*TAU) ' seconds'])
    
    
end
figure;
    NNN=1
    for J=41:5:56
        
        subplot(4,1,NNN)
        plot(X,U(:,J))
        NNN=NNN+1
        set(gca,'YTick',[-1 0 1])
    set(gca,'XTick',[0 0.5 1])
    grid on
        axis([0,1,-1,1])
    %xlabel('X')
    %ylabel('Magnitude')
    %title([['Time = ' num2str((J-1)*TAU)] ' seconds'])
    title([['Time = ' int2str(J-1)] '*tau = ' num2str((J-1)*TAU) ' seconds'])
    
    
    end
figure;
    NNN=1
    for J=61:5:76
        
        subplot(4,1,NNN)
        plot(X,U(:,J))
        set(gca,'YTick',[-1 0 1])
    set(gca,'XTick',[0 0.5 1])
    grid on;
        NNN=NNN+1
        axis([0,1,-1,1])
    %xlabel('X')
    %ylabel('Magnitude')
    %title([['Time = ' num2str((J-1)*TAU)] ' seconds'])
    title([['Time = ' int2str(J-1)] '*tau = ' num2str((J-1)*TAU) ' seconds'])
   
    
end
figure;
    NNN=1
    for J=81:5:96
       
        subplot(4,1,NNN)
        plot(X,U(:,J))
        NNN=NNN+1
        set(gca,'YTick',[-1 0 1])
    set(gca,'XTick',[0 0.5 1])
    grid on;
        axis([0,1,-1,1])
    %xlabel('X')
    %ylabel('Magnitude')
    %title([['Time = ' num2str((J-1)*TAU)] ' seconds'])
    title([['Time = ' int2str(J-1)] '*tau = ' num2str((J-1)*TAU) ' seconds'])
    
    
    end
figure;
    NNN=1
    for J=101:5:116
        
        subplot(4,1,NNN)
        plot(X,U(:,J))
        NNN=NNN+1
        set(gca,'YTick',[-1 0 1])
    set(gca,'XTick',[0 0.5 1])
    grid on;
        axis([0,1,-1,1])
    %xlabel('X')
    %ylabel('Magnitude')
    %title([['Time = ' num2str((J-1)*TAU)] ' seconds'])
    title([['Time = ' int2str(J-1)] '*tau = ' num2str((J-1)*TAU) ' seconds'])
    
    
end