function plotincidentE(E,input)
global nDIM;    

if  nDIM == 2         % Plot 2D contrast/source distribution --------------
    x1 = input.X1(:,1);   x2 = input.X2(1,:);
    set(figure,'Units','centimeters','Position',[5 5 18 12]); 
    subplot(1,2,1);  
        IMAGESC(x1,x2,abs(E{1}));
        title('\fontsize{13} abs(E_1)');  
    subplot(1,2,2);  
        IMAGESC(x1,x2,abs(E{2}));   
        title('\fontsize{13} abs(E_2)');  
      
elseif nDIM == 3      % Plot 3D contrast/source distribution --------------
                      % at x3 = 0 or x3 = dx/2
    N3cross = floor(input.N3/2+1);
    x1 = input.X1(:,1,N3cross);  x2 = input.X2(1,:,N3cross); 
    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,3,1);  
        IMAGESC(x1,x2,abs(E{1}(:,:,N3cross)));
        title('\fontsize{13} abs(E_1)'); 
    subplot(1,3,2);
        IMAGESC(x1,x2,abs(E{2}(:,:,N3cross)));
        title('\fontsize{13} abs(E_2)'); 
end % if 
