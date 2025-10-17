function [data_est] = FunDataSctCircle(P,input)

c_0     = input.c_0;      c_sct   = c_0 / sqrt(1-P(3)); 
gam0    = input.gamma_0;  gam_sct = input.gamma_0 * c_0/c_sct; 

% (1) Compute coefficients of series expansion ----------------------------
  arg0 = gam0 * P(4);  args = gam_sct * P(4); 
  M = 100;                                   % increase M for more accuracy
  A = zeros(1,M+1); 
  for m = 0 : M                     
    Ib0 = besseli(m,arg0);     dIb0 =  besseli(m+1,arg0) + m/arg0 * Ib0; 
    Ibs = besseli(m,args);     dIbs =  besseli(m+1,args) + m/args * Ibs; 
    Kb0 = besselk(m,arg0);     dKb0 = -besselk(m+1,arg0) + m/arg0 * Kb0; 
    A(m+1) = - (gam_sct * dIbs*Ib0 - gam0 * dIb0*Ibs) ...
              /(gam_sct * dIbs*Kb0 - gam0 * dKb0*Ibs);
  end
  
% (2) Compute reflected field at the receivers for all sources (data) -----
  xR(1,:) = input.xR(1,:) - P(1);  % shifted origin 
  xR(2,:) = input.xR(2,:) - P(2);  % shifted origin 
  xS(1,:) = input.xS(1,:) - P(1);  % shifted origin 
  xS(2,:) = input.xS(2,:) - P(2);  % shifted origin 
  rR = sqrt(xR(1,:).^2 + xR(2,:).^2);  phiR = atan2(xR(2,:),xR(1,:));  
  rS = sqrt(xS(1,:).^2 + xS(2,:).^2);  phiS = atan2(xS(2,:),xS(1,:));
    
  data2D = zeros(input.NR,input.NS);
  for p = 1 : input.NR
   for q = 1 : input.NS
    data2D(p,q) = A(1) * besselk(0,gam0*rR(p)) * besselk(0,gam0*rS(q)); 
    for m = 1 : M
     data2D(p,q) = data2D(p,q) + 2 * A(m+1) * cos(m*(phiR(p)-phiS(q))) ...
                    * besselk(m,gam0*rR(p)) * besselk(m,gam0*rS(q)) ;
    end % m_loop   
   end % q_loop
  end % p_loop
  data_est= 1/(2*pi) * data2D;   