function [w_E,dwE] = BICGSTABwdwE(E_inc,dEinc,input)

       CHI = input.CHI_eps;               Rfl = input.Rfl; 
Data_wE{1} = CHI .* E_inc{1};     Data_dwE{1} = Rfl{1} .* dEinc{1};    
Data_wE{2} = CHI .* E_inc{2};     Data_dwE{2} = Rfl{2} .* dEinc{2};

itmax = 200;    it     = 0;     % initialization
w_E{1} = zeros(size(E_inc{1}));        dwE{1} = zeros(size(dEinc{1}));
w_E{2} = zeros(size(E_inc{2}));        dwE{2} = zeros(size(dEinc{2}));
r_E{1} = Data_wE{1};                   drE{1} = Data_dwE{1};
r_E{2} = Data_wE{2};                   drE{2} = Data_dwE{2};             
Norm_D = norm(r_E{1}(:))^2 + norm(r_E{2}(:))^2;
eta_D  = 1 / Norm_D;  Error  = 1;    fprintf('Error =         %g',Error);

while (it < itmax) && ( Error > input.Errcri) %---------------------------% 
  % determine gradient directions 
  AN = sum(r_E{1}(:).*conj(Data_wE{1}(:))+r_E{2}(:).*conj(Data_wE{2}(:)));
  if it == 0;   v_E = r_E;                  dvE = drE;   
  else
     v_E{1} = r_E{1}+(AN/AN_1)*v_E{1};   dvE{1} = drE{1}+(AN/AN_1)*dvE{1};
     v_E{2} = r_E{2}+(AN/AN_1)*v_E{2};   dvE{2} = drE{2}+(AN/AN_1)*dvE{2}; 
  end;                  
  [KvE,KdE] = KopEwdw(v_E,dvE,input); % to determine alpha and beta
  KvE{1} = v_E{1}-CHI.*KvE{1};          KdE{1} = dvE{1}-Rfl{1}.*KdE{1}; 
  KvE{2} = v_E{2}-CHI.*KvE{2};          KdE{2} = dvE{2}-Rfl{2}.*KdE{2}; 
  BN = sum(KvE{1}(:).*conj(Data_wE{1}(:))+KvE{2}(:).*conj(Data_wE{2}(:)));
  alpha = AN / BN;   AN_1   = AN;       
  r_E{1} = r_E{1} - alpha * KvE{1};     drE{1} = drE{1} - alpha * KdE{1};  
  r_E{2} = r_E{2} - alpha * KvE{2};     drE{2} = drE{2} - alpha * KdE{2};
  % + successive overrelaxation (first step of GMRES)
  [KrE,KdrE] = KopEwdw(r_E,drE,input); 
  KrE{1} = r_E{1} - CHI .* KrE{1};   KdrE{1} = drE{1} - Rfl{1} .* KdrE{1};
  KrE{2} = r_E{2} - CHI .* KrE{2};   KdrE{2} = drE{2} - Rfl{2} .* KdrE{2};
  beta = (sum(r_E{1}(:).*conj(KrE{1}(:))+r_E{2}(:).*conj(KrE{2}(:))))... 
         / (norm(KrE{1}(:))^2 +norm(KrE{2}(:))^2);  
  %update contrast sources and residual errors            
  A = alpha; B = beta;  
  w_E{1} = w_E{1}+A*v_E{1}+B*r_E{1};   dwE{1} = dwE{1}+A*dvE{1}+B*drE{1};
  w_E{2} = w_E{2}+A*v_E{2}+B*r_E{2};   dwE{2} = dwE{2}+A*dvE{2}+B*drE{2}; 
  v_E{1} = (A/B)*(v_E{1} - B*KvE{1});  dvE{1} = (A/B)*(dvE{1} - B*KdE{1}); 
  v_E{2} = (A/B)*(v_E{2} - B*KvE{2});  dvE{2} = (A/B)*(dvE{2} - B*KdE{2});
  r_E{1} = r_E{1} - B * KrE{1};        drE{1} = drE{1} - B * KdrE{1};
  r_E{2} = r_E{2} - B * KrE{2};        drE{2} = drE{2} - B * KdrE{2};
  Norm_D = norm(r_E{1}(:))^2  + norm(r_E{1}(:))^2 ; 
  Error=sqrt(eta_D*Norm_D);  fprintf('\b\b\b\b\b\b\b\b%6f',Error); 
   it=it+1;       
end; % CG iterations -----------------------------------------------------
fprintf('\b\b\b\b\b\b\b\b%6f\n',Error); 
disp(['Number of BICGSTABwdwE iterations is ' num2str(it)]);
if it == itmax; disp(['itmax reached: err/norm = ' num2str(Error)]); end;