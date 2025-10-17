function [data_est] = FunDataEMSctCircle(P,input)      
global nDIM;   

% Connect parameters P to the our actual parameters
  eps_sct = P(1);      a = P(2); 

  mu_sct  = input.mu_sct; 
  gam0    = input.gamma_0;         gam_sct = gam0 * sqrt(eps_sct*mu_sct);
  Z_sct   = sqrt(mu_sct/eps_sct); 
 
% (1) Transform Cartesian coordinates to polar coordinates ---------------- 
  xR = input.xR;  
  xS = input.xS;   
  rR = sqrt(xR(1,:).^2 + xR(2,:).^2);  phiR = atan2(xR(2,:),xR(1,:));  
  rS = sqrt(xS(1,:).^2 + xS(2,:).^2);  phiS = atan2(xS(2,:),xS(1,:)); 

% (2) Compute coefficients of Bessel series expansion ---------------------
  arg0 = gam0 * a;  args = gam_sct * a; 
  M = 20;    A = zeros(1,M);                % increase M for more accuracy 
  for m = 1 : M;                      
    Ib0 = besseli(m,arg0);     dIb0 =  besseli(m+1,arg0) + m/arg0 * Ib0; 
    Ibs = besseli(m,args);     dIbs =  besseli(m+1,args) + m/args * Ibs; 
    Kb0 = besselk(m,arg0);     dKb0 = -besselk(m+1,arg0) + m/arg0 * Kb0;
    denominator = Z_sct * dIbs*Kb0 - dKb0*Ibs;
    A(m)  =     -(Z_sct * dIbs*Ib0 - dIb0*Ibs) / denominator; 
  end

% (3) Compute reflected Er field at receivers (data) ----------------------
    data_est = cell(1,nDIM);
  for p = 1 : input.NR;
    for q = 1:input.NS;  
      Er = 0;  Ephi = 0;  ZH3 = 0; 
      for m = 1 : M;
         arg0 = gam0*rR(p);  Kb0 =  besselk(m,arg0);   
                            dKb0 = -besselk(m+1,arg0) + m./arg0 .* Kb0; 
                             KbS =  besselk(m,gam0*rS(q));
           Er = Er  + A(m)*2*m^2.*Kb0.*KbS.*cos(m*(phiR(p)-phiS(q)));
         Ephi = Ephi- A(m)*2*m .*dKb0.*KbS.*sin(m*(phiR(p)-phiS(q)));
         ZH3 = ZH3 - A(m)*2*m  .*Kb0.*KbS.*sin(m*(phiR(p)-phiS(q)));
      end % m_loop
    Er =         1/(2*pi) * Er./rR(p) ./rS(q);
  Ephi =  gam0 * 1/(2*pi) * Ephi      ./rS(q);
  ZH3 = -gam0 * 1/(2*pi) * ZH3       ./rS(q);
    E1 = cos(phiR(p)) .* Er  - sin(phiR(p)) .* Ephi;
    E2 = sin(phiR(p)) .* Er  + cos(phiR(p)) .* Ephi; 
  data_est{1}(p,q)  = (E1); 
  data_est{2}(p,q)  = (E2);
    end % q_loop
  end % p_loop
     
  