function [dataE_est] = FunDataSctIntqwE(P,pG,input)

% Update contrast using analytic contrast profile

          s = 0.01 * input.a;  
          R = sqrt((input.X1).^2 + (input.X2).^2);
    eps_sct = P(1);
 
    input.CHI_eps = (1-eps_sct)./(1+exp((R-P(2))/s));

% Estimate data by solving integral equation with BICGSTABwE method

  dataE_est = cell(1,2);
  for q = 1 : input.NS;  
    [E_incq,~] = IncEMdwave(pG,q,input); 
    [w_E]      = BICGSTABwE(E_incq,input);    % using 3 iterations!
    [Ercv]     = DOPwE(pG,w_E,input); 
    
    dataE_est{1}(:,q) = (Ercv{1});
    dataE_est{2}(:,q) = (Ercv{2});
   end; %q_loop