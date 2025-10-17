function [w_E] = BICGSTABwE(E_inc,input)
 
Data_wE{1}=input.CHI_eps.*E_inc{1};    
Data_wE{2}=input.CHI_eps.*E_inc{2};

itmax = 200;
it     = 0;   % initialization
w_E{1} = zeros(size(E_inc{1}));            w_E{2} = zeros(size(E_inc{2}));
r_E{1} = Data_wE{1};                       r_E{2} = Data_wE{2};
Norm_D = norm(r_E{1}(:))^2 + norm(r_E{2}(:))^2;
eta_D  = 1 / Norm_D;
Error  = 1;            fprintf('Error =         %g',Error); 

while (it < itmax) && ( Error > input.Errcri)
% determine gradient directions 
  AN = sum(r_E{1}(:).*conj(Data_wE{1}(:))+r_E{2}(:).*conj(Data_wE{2}(:)));
     if it == 0 
      v_E{1} = r_E{1};                   v_E{2} = r_E{2};
     else
      v_E{1} = r_E{1}+(AN/AN_1)*v_E{1};  v_E{2} = r_E{2}+(AN/AN_1)*v_E{2}; 
     end    
% determine step length alpha and update residual error      
  KvE     = KopE(v_E,input);
  KvE{1}  = v_E{1} - input.CHI_eps .* KvE{1};
  KvE{2} = v_E{2} - input.CHI_eps .* KvE{2};  
  BN = sum(KvE{1}(:).*conj(Data_wE{1}(:))+KvE{2}(:).*conj(Data_wE{2}(:)));
  alpha = AN / BN;
  r_E{1} = r_E{1} - alpha * KvE{1};  
  r_E{2} = r_E{2} - alpha * KvE{2};
% + succesive overrelaxation (first step of GMRES)
  KrE    = KopE(r_E,input);
  KrE{1} = r_E{1} - input.CHI_eps .* KrE{1};
  KrE{2} = r_E{2} - input.CHI_eps .* KrE{2};
  beta   = (sum(r_E{1}(:).*conj(KrE{1}(:))+r_E{2}(:).*conj(KrE{2}(:))))...
           / (norm(KrE{1}(:))^2 + norm(KrE{2}(:))^2);
% update contrast sources and w , v and AN 
  w_E{1} = w_E{1} + alpha * v_E{1} + beta * r_E{1}; 
  w_E{2} = w_E{2} + alpha * v_E{2} + beta * r_E{2}; 
  v_E{1} = (alpha/beta)  * (v_E{1} - beta * KvE{1});
  v_E{2} = (alpha/beta)  * (v_E{2} - beta * KvE{2});
  AN_1   = AN;
% update residual errors r_D  
  r_E{1} = r_E{1} - beta * KrE{1};
  r_E{2} = r_E{2} - beta * KrE{2};
  Norm_D = norm(r_E{1}(:))^2 + norm(r_E{2}(:))^2; 
  Error  = sqrt(eta_D * Norm_D);    % fprintf('\b\b\b\b\b\b\b\b%6f',Error);        
   it = it+1;
 end % CG iterations
 fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
 disp(['Number of BICGSTABwE iterations is ' num2str(it)]);
 if it == itmax; disp(['itmax reached: err/norm = ' num2str(Error)]); end