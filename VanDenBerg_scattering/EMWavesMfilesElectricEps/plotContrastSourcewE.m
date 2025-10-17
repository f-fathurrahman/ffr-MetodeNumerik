function plotContrastSourcewE(w_E,input)
global nDIM;     set(figure,'Units','centimeters','Position',[5 5 18 12]);
 
if  nDIM == 2       % Plot 2D contrast/source distribution ---------------
    x1 = input.X1(:,1);   x2 = input.X2(1,:);
    subplot(1,2,1);  
        IMAGESC(x1,x2,abs(w_E{1}));
        title('\fontsize{13} abs(w_1^E)');  
    subplot(1,2,2);  
        IMAGESC(x1,x2,abs(w_E{2}));   
        title('\fontsize{13} abs(w_2^E)');          
elseif nDIM == 3    % Plot 3D contrast/source distribution ---------------
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1); 
    N3cross = floor(input.N3/2+1);           
    subplot(1,3,1);
        IMAGESC(x1,x2,abs(w_E{1}(:,:,N3cross)));
        title('\fontsize{13} abs(w_1^E)');   
    subplot(1,3,2);
        IMAGESC(x1,x2,abs(w_E{2}(:,:,N3cross)));
        title('\fontsize{13} abs(w_2^E)');  
    subplot(1,3,3);
       IMAGESC(x1,x2,abs(w_E{3}(:,:,N3cross)));
        title('\fontsize{13} abs(w_3^E)'); 
end % if 