function [w_E,dwE] = BiCGSTABwdwE(E_inc,dEinc,input)

 CHI = input.CHI_eps;  
 Rfl = input.Rfl; 
Data_wE{1}  = CHI .* E_inc{1};    
Data_wE{2}  = CHI .* E_inc{2};
Data_dwE{1} = Rfl{1} .* dEinc{1};    
Data_dwE{2} = Rfl{2} .* dEinc{2};

itmax = 200;
it     = 0;   % initialization
w_E{1} = zeros(size(E_inc{1}));            
w_E{2} = zeros(size(E_inc{2}));
dwE{1} = zeros(size(dEinc{1}));            
dwE{2} = zeros(size(dEinc{2}));
r_E{1} = Data_wE{1};                       
r_E{2} = Data_wE{2};
drE{1} = Data_dwE{1};                       
drE{2} = Data_dwE{2};
Norm_D =  norm(r_E{1}(:))^2 + norm(r_E{2}(:))^2  ...
        + norm(drE{1}(:))^2 + norm(drE{2}(:))^2;
eta_D  = 1 / Norm_D;
Error  = 1;            fprintf('Error =         %g',Error); 

while (it < itmax) && ( Error > input.Errcri)
% determine gradient directions 
  AN =  sum(r_E{1}(:).*conj(Data_wE{1}(:))) ...
      + sum(r_E{2}(:).*conj(Data_wE{2}(:))) ...
      + sum(drE{1}(:).*conj(Data_dwE{1}(:))) ...
      + sum(drE{2}(:).*conj(Data_dwE{2}(:)));
     if it == 0; 
      v_E{1} = r_E{1};                   
      v_E{2} = r_E{2};
      dvE{1} = drE{1};                   
      dvE{2} = drE{2};
     else
      v_E{1} = r_E{1} + (AN/AN_1) * v_E{1};
      v_E{2} = r_E{2} + (AN/AN_1) * v_E{2}; 
      dvE{1} = drE{1} + (AN/AN_1) * dvE{1};
      dvE{2} = drE{2} + (AN/AN_1) * dvE{2};
     end    
% determine step length alpha and update residual error      
  [KvE,KdvE] = KopEwdw(v_E,dvE,input);
  KvE{1}  = v_E{1} - CHI .* KvE{1};
  KvE{2}  = v_E{2} - CHI .* KvE{2};  
  KdvE{1} = dvE{1} - Rfl{1} .* KdvE{1};
  KdvE{2} = dvE{2} - Rfl{2} .* KdvE{2}; 
  
  BN =  sum(KvE{1}(:) .* conj(Data_wE{1}(:)))   ...
      + sum(KvE{2}(:) .* conj(Data_wE{2}(:)))   ... 
      + sum(KdvE{1}(:) .* conj(Data_dwE{1}(:))) ...
      + sum(KdvE{2}(:) .* conj(Data_dwE{2}(:)));
  alpha = AN / BN;
  r_E{1} = r_E{1} - alpha * KvE{1};  
  r_E{2} = r_E{2} - alpha * KvE{2};
  drE{1} = drE{1} - alpha * KdvE{1};  
  drE{2} = drE{2} - alpha * KdvE{2};
  
% + succesive overrelaxation (first step of GMRES)
  [KrE,KdrE] = KopEwdw(r_E,drE,input);
  KrE{1}  = r_E{1} - CHI .* KrE{1};
  KrE{2}  = r_E{2} - CHI .* KrE{2};
  KdrE{1} = drE{1} - Rfl{1} .* KdrE{1};
  KdrE{2} = drE{2} - Rfl{2} .* KdrE{2};
  beta    = (  sum(r_E{1}(:).*conj(KrE{1}(:)))    ...
             + sum(r_E{2}(:).*conj(KrE{2}(:)))    ...
             + sum(drE{1}(:).*conj(KdrE{1}(:)))   ...
             + sum(drE{2}(:).*conj(KdrE{2}(:))) ) ... 
             / (  norm(KrE{1}(:))^2  + norm(KrE{2}(:))^2 ...
                + norm(KdrE{1}(:))^2 + norm(KdrE{2}(:))^2 );
% update contrast sources and w , v and AN 
  w_E{1} = w_E{1} + alpha * v_E{1} + beta * r_E{1}; 
  w_E{2} = w_E{2} + alpha * v_E{2} + beta * r_E{2}; 
  dwE{1} = dwE{1} + alpha * dvE{1} + beta * drE{1}; 
  dwE{2} = dwE{2} + alpha * dvE{2} + beta * drE{2}; 
  v_E{1} = (alpha/beta)  * (v_E{1} - beta * KvE{1});
  v_E{2} = (alpha/beta)  * (v_E{2} - beta * KvE{2});
  dvE{1} = (alpha/beta)  * (dvE{1} - beta * KdvE{1});
  dvE{2} = (alpha/beta)  * (dvE{2} - beta * KdvE{2});
  AN_1   = AN;
% update residual errors r_D  
  r_E{1} = r_E{1} - beta * KrE{1};
  r_E{2} = r_E{2} - beta * KrE{2};
  drE{1} = drE{1} - beta * KdrE{1};
  drE{2} = drE{2} - beta * KdrE{2};
  Norm_D =  norm(r_E{1}(:))^2 + norm(r_E{2}(:))^2 ...
          + norm(drE{1}(:))^2 + norm(drE{2}(:))^2 ; 
  Error  = sqrt(eta_D * Norm_D);     fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
   it = it+1;
 end; % CG iterations
 fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
 disp(['Number of iterations is ' num2str(it)]);
 if it == itmax; disp(['itmax reached: err/norm = ' num2str(Error)]); end;