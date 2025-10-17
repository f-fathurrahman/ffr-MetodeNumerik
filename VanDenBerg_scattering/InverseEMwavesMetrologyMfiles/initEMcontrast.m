function [input] = initEMcontrast(input);
input.a  = 40;      % radius circle cylinder 
       R = sqrt(input.X1.^2 + input.X2.^2);

input.xO(1) = -input.a * 0.4;   % center coordinates of two small circles
input.xO(2) = input.a  * 0.5 ; 
         R1 = sqrt((input.X1-input.xO(1)).^2 + (input.X2-input.xO(2)).^2);        
input.xO(1) = -input.a * 0.4;   % center coordinates of circle
input.xO(2) = -input.a * 0.5; 
         R2 = sqrt((input.X1-input.xO(1)).^2 + (input.X2-input.xO(2)).^2); 

% Make a choice of contrast profiles for the base model
  WindowCircle= (R < input.a); 

  Window_eyes = (R1 < 0.1*input.a) +  (R2 < 0.1*input.a);

  Window_nose = (abs(input.X1)<0.2*input.a).*(abs(input.X2)<0.08*input.a); 

  Window_mouth= (abs(input.X1-0.5*input.a) < 0.05*input.a) ...
                                          .*(abs(input.X2) < 0.5*input.a);               
    
 % Set factors either to one or zero, for presence/absence of model                                      
   A_eyes  = 1;     A_nose  = 1;   A_mouth = 1;     
   input.CHI_eps = (1-input.eps_sct)  * (WindowCircle            ...
                                      - A_eyes  * Window_eyes/2  ...
                                      - A_nose  * Window_nose/2  ...
                                      - A_mouth * Window_mouth/2);
% Plot base profile
  figure(1); 
  x1 = input.X1(:,1); x2 = input.X2(1,:); 
  IMAGESC(x1,x2,real(input.CHI_eps));    colormap(hot); grid; grid minor; 