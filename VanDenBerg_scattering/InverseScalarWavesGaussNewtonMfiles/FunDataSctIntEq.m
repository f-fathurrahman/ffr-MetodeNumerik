function [data_est] = FunDataSctIntEq(P,G_S,G_R,input)
w = cell(1,input.NS); data_est = zeros(input.NR,input.NS);

% Update contrast using analytic contrast profile
          s = 0.3 * input.dx;
          R = sqrt((input.X1-P(1)).^2 + (input.X2-P(2)).^2);
  input.CHI = P(3)./(1+exp((R-P(4))/s));

% Solve integral equation for contrast source with BICGSTAB method -------            
  for q = 1 : input.NS             
      w{q} = BICGSTABw(G_S{q},input);    
  end % q_loop

% Compute synthetic data and plot fields and data ------------------------
  for q = 1 : input.NS
      data_est(:,q) = DopM(G_R,w{q},input);
  end % q_loop