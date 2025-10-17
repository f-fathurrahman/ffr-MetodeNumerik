function [dataE_est] = FunDataSctIntqwdwE(P,pG,input)

% Update contrast using analytic contrast and determine interface contrast
              s = 0.01 * input.a;  
              R = sqrt(input.X1.^2 + input.X2.^2);      [N1,N2] = size(R);
  input.CHI_eps = (1-P(1))./(1+exp((R-P(2))/s));  
            eps = 1 - input.CHI_eps; 
            
          Rfl{1} = zeros(N1,N2);      Rfl{2} = zeros(N1,N2); 
Rfl{1}(1:N1-1,:) = (eps(2:N1,:) - eps(1:N1-1,:)) ...
                                      ./(eps(2:N1,:) + eps(1:N1-1,:));   
Rfl{2}(:,1:N2-1) = (eps(:,2:N2) - eps(:,1:N2-1)) ...
                                      ./(eps(:,2:N2) + eps(:,1:N2-1));
input.Rfl = Rfl;  

% Estimate data by solving integral equation with BICGSTABwE method 
dataE_est = cell(1,2);       
  for q = 1 : input.NS;  
      [E_incq,dEincq] = IncEMdwave(pG,q,input); 
      [w_E,dwE] = BICGSTABwdwE(E_incq,dEincq,input); % using 3 iterations!
      [Ercv]    = DOPEwdw(w_E,dwE,input); 
      dataE_est{1}(:,q) = (Ercv{1});
      dataE_est{2}(:,q) = (Ercv{2});
  end; %q_loop