function [Pinc] = IncPressureWave(q,input) 
gam0 = input.gamma_0;  
xS   = input.xS;    dx = input.dx;
X1     = input.X1;   X2 = input.X2; 

DIS     = sqrt( (X1+dx/2-xS(1,q)).^2 + (X2-xS(2,q)).^2 );
Pinc{1} = 1/(2*pi).* besselk(0,gam0*DIS);   % Grid (1) on shifted points               

DIS     = sqrt( (X1-xS(1,q)).^2 + (X2+dx/2-xS(2,q)).^2 );
Pinc{2} = 1/(2*pi).* besselk(0,gam0*DIS);   % Grid (2) on shifted points  