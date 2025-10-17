function [w,u,CHI,r_M,r_D,eta_D,ErrorD,ErrorM] = ...
                         InitialEstimate(eta_M,u_inc,G_R,data,input)
global it                    
FFTG = input.FFTG;  NS  = input.NS;   
 w   = cell(1,NS);   u  = cell(1,NS);   r_M = cell(1,NS);  

% Determine contrast sources by back-projection and minimization of Norm_M
  it = 0;
  Norm_M = 0; 
  for q = 1 : NS
    r_M{q} = data(:,q);
    g      = AdjDopM(G_R,r_M{q},input);                
    A      = norm(g(:))^2;
    GRv    = DopM(G_R,g,input);               
    B      = norm(GRv(:))^2; 
    alpha  = real(A/B);
    r_M{q} = r_M{q} - alpha * GRv(:);
    w{q}   = alpha * g; 
    u{q}   = u_inc{q} + Kop(w{q},FFTG);
   Norm_M  = Norm_M + norm(r_M{q}(:))^2; 
  end %q_loop
  ErrorM = eta_M * Norm_M;

% Update contrast  by  minimization of  Norm_D = || CHI * u_inc - w ||^2 
  [CHI,r_D,eta_D,ErrorD] = UpdateContrast(w,u_inc,u_inc,ErrorM,input);
  
  disp([' ErrorM = ',num2str(ErrorM),'   ErrorD = ',num2str(ErrorD)]);
  disp(['Iteration = ', num2str(it)]);                      saveCHI(CHI);
  disp('--------------------------------------------------------------');

% Figures of contrast
  x1 = input.X1(:,1);   x2 = input.X2(1,:); 
  figure(3);  IMAGESC(x1,x2,real(CHI));  
  figure(4);  mesh(abs(CHI),'edgecolor', 'k'); 
              view(37.5,45); axis xy; axis('off'); axis('tight') 
              axis([1 input.N1 1 input.N2  0 1.5*max(abs(input.CHI(:)))]);
end