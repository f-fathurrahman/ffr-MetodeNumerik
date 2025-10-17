function [Rfl,dr,eta_D,ErrorD] ...
                    = UpdateContrastdw(dw,P,P_inc,ErrorM,input)
global nDIM
dr  = cell(nDIM,input.NS);          
% Determine mass density interface contrast
  INP1_dw1 = zeros(input.N1,input.N2);     
  INP1_P1  = zeros(input.N1,input.N2);  
  for q = 1 : input.NS     
    INP1_dw1 = INP1_dw1 + conj(P{1,q}) .* dw{1,q};
    INP1_P1  = INP1_P1  + conj(P{1,q}) .*  P{1,q};
  end %q_loop 
  Rfl{1} = INP1_dw1 ./ INP1_P1; 
   shape = abs(Rfl{1}) > 0.10 * max(abs(Rfl{1}(:)));   
  Rfl{1} = shape .* Rfl{1};              % Neglect small small variations
 
  INP2_dw2 = zeros(input.N1,input.N2);     
  INP2_P2  = zeros(input.N1,input.N2);  
  for q = 1 : input.NS     
    INP2_dw2 = INP2_dw2 + conj(P{2,q}) .* dw{2,q};
    INP2_P2  = INP2_P2  + conj(P{2,q}) .* P{2,q};
  end %q_loop 
  Rfl{2} = INP2_dw2 ./ INP2_P2;
   shape = abs(Rfl{2}) >  0.10 * max(abs(Rfl{2}(:)));    
  Rfl{2} = shape .* Rfl{2};              % Neglect small small variations
  
% Compute residual error in object domain 
  Norm_D = 0;   Norm = 0;  
  for q = 1 : input.NS
    dr{1,q}  = Rfl{1} .*  P{1,q}   - dw{1,q};  
    dr{2,q}  = Rfl{2} .*  P{2,q}   - dw{2,q};   
    Norm_D = Norm_D + norm(dr{1,q}(:))^2 +  norm(dr{2,q}(:))^2 ;
    Norm   = Norm   + norm(Rfl{1}(:).*P_inc{1,q}(:))^2 ...
                    + norm(Rfl{2}(:).*P_inc{2,q}(:))^2;
  end %q_loop 
  eta_D = 1 / Norm;
  ErrorD = eta_D * Norm_D;
  disp(['  ErrorM = ',num2str(ErrorM), '    ErrorD  = ',num2str(ErrorD)]);

% Compute the mass-density volume contrast from interface contrast
  CHI_rho = RfltoCHI(Rfl,input);  

% Show intermediate pictures of reconstruction of compressibility contrast
  x1 = input.X1(:,1);   x2 = input.X2(1,:);   
  set(figure(6),'Units','centimeters','Position',[5 1 27 12]); 
  subplot(1,3,1); IMAGESC(x1+input.dx/2,x2,real(Rfl{1}));
                  title('\fontsize{13} R^{(1)}');
  subplot(1,3,2); IMAGESC(x1,x2+input.dx/2,real(Rfl{2}));
                  title('\fontsize{13} R^{(2)}');   
  subplot(1,3,3); IMAGESC(x1,x2,real(CHI_rho));  
                  title('\fontsize{13} 1-\rho_{sct}/ \rho_{0}'); pause(.1)  