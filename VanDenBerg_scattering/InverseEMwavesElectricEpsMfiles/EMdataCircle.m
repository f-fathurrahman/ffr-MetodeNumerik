  input = initEM();       
      a = input.a;                       gam0    = input.gamma_0; 
eps_sct = input.eps_sct;                 mu_sct  = input.mu_sct;
gam_sct = gam0 * sqrt(eps_sct*mu_sct);   Z_sct   = sqrt(mu_sct/eps_sct);
    
% (1) Transform Cartesian coordinates to polar coordinates ---------------- 
  xR = input.xR;  
  xS = input.xS;   
  rR = sqrt(xR(1,:).^2 + xR(2,:).^2);  phiR = atan2(xR(2,:),xR(1,:));  
  rS = sqrt(xS(1,:).^2 + xS(2,:).^2);  phiS = atan2(xS(2,:),xS(1,:)); 

% (2) Compute coefficients of Bessel series expansion ---------------------
  arg0 = gam0 * input.a;  args = gam_sct * input.a; 
  M = 20;    A = zeros(1,M);                % increase M for more accuracy 
  for m = 1 : M                      
    Ib0 = besseli(m,arg0);     dIb0 =  besseli(m+1,arg0) + m/arg0 * Ib0; 
    Ibs = besseli(m,args);     dIbs =  besseli(m+1,args) + m/args * Ibs; 
    Kb0 = besselk(m,arg0);     dKb0 = -besselk(m+1,arg0) + m/arg0 * Kb0;
    denominator = Z_sct * dIbs*Kb0 - dKb0*Ibs;
    A(m)  =     -(Z_sct * dIbs*Ib0 - dIb0*Ibs) / denominator; 
  end
  
% (3) Compute reflected Er field at receivers (data) ----------------------
  E1data = zeros(input.NR,input.NS); E2data =  zeros(input.NR,input.NS);
  for p = 1 : input.NR
    for q = 1:input.NS  
      Er = 0;  Ephi = 0; 
      for m = 1 : M
         arg0 = gam0*rR(p);  Kb0 =  besselk(m,arg0);   
                            dKb0 = -besselk(m+1,arg0) + m./arg0 .* Kb0; 
                             KbS =  besselk(m,gam0*rS(q));
           Er = Er  + A(m)*2*m^2.*Kb0.*KbS.*cos(m*(phiR(p)-phiS(q)));
         Ephi = Ephi- A(m)*2*m .*dKb0.*KbS.*sin(m*(phiR(p)-phiS(q)));
      end % m_loop
    Er =         1/(2*pi) * Er./rR(p) ./rS(q);
  Ephi =  gam0 * 1/(2*pi) * Ephi      ./rS(q);
    E1 = cos(phiR(p)) .* Er  - sin(phiR(p)) .* Ephi;
    E2 = sin(phiR(p)) .* Er  + cos(phiR(p)) .* Ephi;
  E1data(p,q) = E1;  E2data(p,q) = E2;   
    end % q_loop
  end % p_loop                                 
  save Edata E1data E2data;        