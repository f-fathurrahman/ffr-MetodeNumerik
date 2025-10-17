clear all; clc; close all; clear workspace
global Pnom P_it
input = initEM();  itmaxGN = 15;     eps = 1e-6;

% Compute Green functions for data operators 
  pG = EGreenSourceReceivers(input);

% Computes analytical data for actual parameters
  P = [input.eps_sct, input.a];  EMdata = FunDataEMSctCircle(P,input);

% Start with estimated parameters
  Pnom = [input.eps_sct,input.a];  disp(['Nominal P : ', num2str(Pnom)]);
  P0   = [1.05,         70     ];  disp(['Initial P0: ', num2str(P0)]);
  P_it(1,:) = P0;   % save as global array for plotting
 
 for it = 1: itmaxGN 
  % Compute residual error, construct Jacobian matrix and update P
    P = P_it(it,:);  
    data_est    = FunDataSctIntqwE(P,pG,input);
    Residual{1} = EMdata{1} -  data_est{1};  
    Residual{2} = EMdata{2} -  data_est{2}; 
    GNerror     = (norm(Residual{1}(:))+norm(Residual{2}(:)))...
                 /(norm(EMdata{1}(:))+norm(EMdata{2}(:))) ;
    disp(['           it= ',num2str(it),'  Residual= ',num2str(GNerror)]);
  % Compute numerical derivatives with respect to P   
    P(1) = P(1) + eps;  data_plus = FunDataSctIntqwE(P,pG,input);
      dCHI1_P{1} = (data_plus{1} - data_est{1}) /eps;  
      dCHI2_P{1} = (data_plus{2} - data_est{2}) /eps;  
    P(1) = P(1) - eps;
    P(2) = P(2) + eps;  data_plus = FunDataSctIntqwE(P,pG,input);
     dCHI1_P{2} = (data_plus{1} - data_est{1}) /eps;
     dCHI2_P{2} = (data_plus{2} - data_est{2}) /eps;       
    P(2) = P(2) - eps;
  % Compute residual error, construct Jacobian matrix and update P
    b = zeros(1,2); A = zeros(2,2);
    for l = 1:2;
        b(l) = sum(conj(dCHI1_P{l}(:)).*Residual{1}(:))...
             + sum(conj(dCHI2_P{l}(:)).*Residual{2}(:));
        for k = 1:2;
            A(l,k) = sum(conj(dCHI1_P{l}(:)).*dCHI1_P{k}(:))...
                   + sum(conj(dCHI2_P{l}(:)).*dCHI2_P{k}(:));
        end;
    end;                                 
    DeltaP = real(b/A);           alpha = 1;   
    alpha  = 0.618;            %  use Golden ratio for damping  
    P = P + alpha * DeltaP ;         
    disp(['        P: ', num2str(P)]);
    P_it(it+1,:) = P; % save as global array for plotting
 end; % it_loop
 
 DisplayParameterswE;