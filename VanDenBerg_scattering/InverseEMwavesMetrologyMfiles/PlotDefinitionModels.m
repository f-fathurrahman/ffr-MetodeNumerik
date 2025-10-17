input = initEM();
x1 = input.X1(:,1);   x2 = input.X2(1,:);   
input.a  = 40;      % radius circle cylinder 
       R = sqrt(input.X1.^2 + input.X2.^2);

input.xO(1) = -input.a * 0.4;   % center coordinates of circle
input.xO(2) = input.a  * 0.5 ; 
          R1 = sqrt((input.X1-input.xO(1)).^2 + (input.X2-input.xO(2)).^2);        
input.xO(1) = -input.a * 0.4;   % center coordinates of circle
input.xO(2) = -input.a * 0.5; 
          R2 = sqrt((input.X1-input.xO(1)).^2 + (input.X2-input.xO(2)).^2);
          
% Make contrast profiles with various defects
window_Circle = (R < input.a); 

window_Eyes   = (R1 < 0.1*input.a) +  (R2 < 0.1*input.a);

window_Nose   = (abs(input.X1)<0.2*input.a).*(abs(input.X2)<0.08*input.a); 

window_Mouth  = (abs(input.X1-0.5*input.a) < 0.05*input.a) ...
                                          .*(abs(input.X2) < 0.5*input.a);               
    
% Set factors either to one or zero, for presence/absence of model 
CHI_eps1 = (1-input.eps_sct)  * window_Circle;

A_Eyes  = 1;     A_Nose  = 1;   A_Mouth = 1;     
CHI_eps2 = (1-input.eps_sct)  * (window_Circle              ...
                                      - A_Eyes  * window_Eyes/2  ...
                                      - A_Nose  * window_Nose/2  ...
                                      - A_Mouth * window_Mouth/2);        
               
figure; 
subplot(1,3,1);
  IMAGESC(x1,x2,real(CHI_eps2));                 
  title('\fontsize{11}Base model'); colorbar off; axis off;
   xlabel(''); ylabel('');  
subplot(1,3,2);
  IMAGESC(x1,x2,real(CHI_eps1));                 
  title('\fontsize{11}Monitor model');  colorbar off;axis off;
    xlabel(''); ylabel('');     
subplot(1,3,3)
  IMAGESC(x1,x2,real(CHI_eps1-CHI_eps2));                 
  title('\fontsize{11}Defects model');  colorbar off;axis off;
   xlabel(''); ylabel(''); 
 colormap(jet.^(1/4));