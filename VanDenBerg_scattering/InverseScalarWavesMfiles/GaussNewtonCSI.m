function [CHI,input] = GaussNewtonCSI(CHI_csi,input)
global it Pnom P_it
X1 = input.X1;  X2 = input.X2;

% Start with estimated contrast sources
if it == 1 
   Pnom = [input.xO(1), input.xO(2), input.contrast,   input.a];  
   P0   = [0, 0, 2*mean(CHI_csi(:)), 75];
   disp(['Nominal P: ', num2str(Pnom)]);
   disp(['Initial P: ', num2str(P0)]);
   P_it(1,:) = P0;   % save as global array for plotting
end
  
  P = P_it(it,:);  
% Determine weak function and for a typical value of s  
  s = 2 * input.dx;  if it > 100;  s =  1 * input.dx;  end  
                     
  R = sqrt((X1-P(1)).^2 + (X2-P(2)).^2);       Rp = R - P(4) + input.dx/2;
                                               Rm = R - P(4) - input.dx/2;
  CHI_weak = -P(3) * (log(1+exp(-Rp/s))-log(1+exp(-Rm/s))) / (input.dx/s);
  dCHI_R   =  P(3) * (1./(1+exp(Rp/s)) - 1./(1+exp(Rm/s))) /  input.dx;
  
  if sum(isnan(CHI_weak(:))) > 0 
       disp(['parameter s too small: increase value of ',num2str(s)]);
  end
  
 % Compute analytical derivatives with respect to P   
  dCHI_P{1} = -(X1-P(1))./R .* dCHI_R;
  dCHI_P{2} = -(X2-P(2))./R .* dCHI_R;
  dCHI_P{3} = CHI_weak / P(3);
  dCHI_P{4} = -dCHI_R;
  
% Compute residual error, construct Jacobian matrix and update P
  Residual = CHI_csi - CHI_weak;  
  GNerror  = norm(Residual(:))/norm(CHI_csi(:));
  disp(['           it= ',num2str(it),'  Residual= ',num2str(GNerror)]);
  b = zeros(1,4); A = zeros(4,4);
  for l = 1:4
      b(l) = sum(dCHI_P{l}(:).*Residual(:));
      for k = 1:4
          A(l,k) = sum(dCHI_P{l}(:).*dCHI_P{k}(:));
      end
  end                                 
  DeltaP = real(b/A);  P = P + DeltaP; P(4) = abs(P(4)); % P(4)is positive 
  disp(['        P: ', num2str(P)]);
  P_it(it+1,:) = P;  % save as global array for plotting
  
% Update contrast using analytic contrast profile
  R   = sqrt((X1-P(1)).^2 + (X2-P(2)).^2);
  CHI = P(3)./(1+exp((R-P(4))/s));