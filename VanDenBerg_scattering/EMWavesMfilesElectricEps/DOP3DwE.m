function [Edata,Hdata] = DOP3DwE(w_E,input)
gam0 = input.gamma_0;    dx = input.dx;   xR = input.xR; 

Edata = zeros(1,input.NR); 
Hdata = zeros(1,input.NR); 

% Weak form
delta  = (4*pi/3)^(-1/3) * dx;       % radius sphere with area of dx^3
   arg = gam0*delta; 
factor = 3 * (cosh(arg) - sinh(arg)/arg) / arg^2; 

for p = 1 : input.NR
    X1  = xR(1,p)-input.X1;  X2 = xR(2,p)-input.X2; X3 = xR(3,p)-input.X3;
    DIS = sqrt(X1.^2+X2.^2+X3.^2); X1=X1./DIS;  X2=X2./DIS;  X3=X3./DIS;
    G   = factor * exp(-gam0*DIS) ./ (4*pi*DIS);    

    dG  =  - (gam0 + 1./DIS) .* G;  
   d1_G = X1 .* dG;       d2_G = X2 .* dG;       d3_G = X3 .* dG; 
   
   dG11 = ( (3*X1.*X1 - 1).*(gam0./DIS+1./DIS.^2) + gam0^2*X1.*X1 ) .* G;
   dG22 = ( (3*X2.*X2 - 1).*(gam0./DIS+1./DIS.^2) + gam0^2*X2.*X2 ) .* G;
   dG33 = ( (3*X3.*X3 - 1).*(gam0./DIS+1./DIS.^2) + gam0^2*X3.*X3 ) .* G;
   dG21 = (  3*X2.*X1     .*(gam0./DIS+1./DIS.^2) + gam0^2*X2.*X1 ) .* G;
   dG31 = (  3*X3.*X1     .*(gam0./DIS+1./DIS.^2) + gam0^2*X3.*X1 ) .* G;
   dG32 = (  3*X3.*X2     .*(gam0./DIS+1./DIS.^2) + gam0^2*X3.*X2 ) .* G;
  
  E1rfl = dx^3 * sum( (gam0^2*G(:) - dG11(:)) .* w_E{1}(:)  ...
                                   - dG21(:)  .* w_E{2}(:)  ...
                                   - dG31(:)  .* w_E{3}(:) ); 
  E2rfl = dx^3 * sum(              - dG21(:)  .* w_E{1}(:)  ...
                     +(gam0^2*G(:) - dG22(:)) .* w_E{2}(:)  ...
                                   - dG32(:)  .* w_E{3}(:) );  
  E3rfl = dx^3 * sum(              - dG31(:)  .* w_E{1}(:)  ...
                                   - dG32(:)  .* w_E{2}(:)  ...
                     +(gam0^2*G(:) - dG33(:)) .* w_E{3}(:) ); 
 
 ZH1rfl = -gam0 * dx^3 * sum( d2_G(:).*w_E{3}(:) - d3_G(:).*w_E{2}(:) );
 ZH2rfl = -gam0 * dx^3 * sum( d3_G(:).*w_E{1}(:) - d1_G(:).*w_E{3}(:) );
 ZH3rfl = -gam0 * dx^3 * sum( d1_G(:).*w_E{2}(:) - d2_G(:).*w_E{1}(:) );
 
 Edata(1,p) = sqrt(abs(E1rfl)^2  + abs(E2rfl)^2  + abs(E3rfl^2));
 Hdata(1,p) = sqrt(abs(ZH1rfl)^2 + abs(ZH2rfl)^2 + abs(ZH3rfl^2));
 
end  % p_loop