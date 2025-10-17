function [DivFocus] = divFocus(data1,data2,input)
  N1 = input.N1; N2 = input.N2;
  
% Compute receiver focusing operators ------------------------------------
  Focus_R = cell(1,input.NR); 
  for p = 1 : input.NR     
        X1  = input.X1 - input.xR(1,p); 
        X2  = input.X2 - input.xR(2,p);   
        DIS = sqrt(X1.^2 + X2.^2);  
        Focus_R{p} = 1./ besselk(0,input.gamma_0*DIS);
  end; % p_loop

% Compute source focusing operators --------------------------------------
  Focus_S = cell(1,input.NS);
  for q = 1 : input.NS     
        X1 = input.X1 - input.xS(1,q); 
        X2 = input.X2 - input.xS(2,q);   
       DIS = sqrt(X1.^2 + X2.^2);  
       Focus_S{q} = 1./ besselk(0,input.gamma_0*DIS); 
  end; % q_loop

% Focus the two components of the electric field data --------------------
  dataFocus{1} = zeros(input.N1,input.N2);
  dataFocus{2} = zeros(input.N1,input.N2); 
  for p = 1:input.NR;
    for q = 1 :input.NS;
     dataFocus{1} = dataFocus{1} + data1(p,q) .* Focus_R{p} .* Focus_S{q}; 
     dataFocus{2} = dataFocus{2} + data2(p,q) .* Focus_R{p} .* Focus_S{q};
    end; % q_loop  
  end; % p_loop

% Take divergence of vector (dataFocus{1},dataFocus{2}) ------------------
  DivFocus = DIV(dataFocus,input);
  % Enforce the divergence at the boundary of the test domain to zero
    DivFocus(1,:)  = DivFocus(2,:);    DivFocus(:,1)  = DivFocus(:,2); 
    DivFocus(N1,:) = DivFocus(N1-1,:); DivFocus(:,N2) = DivFocus(:,N2-1);