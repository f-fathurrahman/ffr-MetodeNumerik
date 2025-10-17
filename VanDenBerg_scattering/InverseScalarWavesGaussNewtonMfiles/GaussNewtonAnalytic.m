clear all; clc; close all; clear workspace
global Pnom P_it
input = init();
itmax = 15;
eps = 1e-6;

% (1) Compute Exact data using analytic Bessel function expansion -------
      disp('Running DataSctCircle');  DataSctCircle;  
      load DataAnalytic;
        
% (2) Start with estimated contrast sources
   Pnom = [input.xO(1), input.xO(2), input.contrast,   input.a];  
   P0   = [ 0, 0, input.contrast/10, 75];
   disp(['Nominal P : ', num2str(Pnom)]);
   disp(['Initial P0: ', num2str(P0)]);
   P_it(1,:) = P0;   % save as global array for plotting
 
 for it = 1: itmax 
  % Compute residual error, construct Jacobian matrix and update P
    P = P_it(it,:);  
    data_est = FunDataSctCircle(P,input);
    Residual = dataCircle -  data_est;  
    GNerror  = norm(Residual(:))/norm(dataCircle(:));
    disp(['           it= ',num2str(it),'  Residual= ',num2str(GNerror)]);
  
  % Compute numerical derivatives with respect to P   
    P(1) = P(1) + eps;  data_plus = FunDataSctCircle(P,input);
      dCHI_P{1} = (data_plus - data_est) /eps;       P(1) = P(1) - eps;
    P(2) = P(2) + eps;  data_plus = FunDataSctCircle(P,input);
      dCHI_P{2} = (data_plus - data_est) /eps;       P(2) = P(2) - eps;
    P(3) = P(3) + eps;  data_plus = FunDataSctCircle(P,input);
      dCHI_P{3} = (data_plus - data_est) /eps;       P(3) = P(3) - eps; 
    P(4) = P(4) + eps;  data_plus = FunDataSctCircle(P,input);
      dCHI_P{4} = (data_plus - data_est) /eps;       P(4) = P(4) - eps;  
   
  % Compute residual error, construct Jacobian matrix and update P
    b = zeros(1,4); A = zeros(4,4);
    for l = 1:4
        b(l) = sum(conj(dCHI_P{l}(:)).*Residual(:));
        for k = 1:4
            A(l,k) = sum(conj(dCHI_P{l}(:)).*dCHI_P{k}(:));
        end
    end                                 
    DeltaP = real(b/A); 
    P = P + DeltaP;              P(4) = abs(P(4)); % P(4) is positive 
    disp(['        P: ', num2str(P)]);
    P_it(it+1,:) = P; % save as global array for plotting
 end % it_loop
 
 DisplayParameters;