function plotEMInterfaceSourceE(dwE,input)
global nDIM;  
dx = input.dx;       Rfl = input.Rfl;

if  nDIM == 2
    x1 = input.X1(:,1);   x2 = input.X2(1,:);
    
    set(figure,'Units','centimeters','Position',[5 5 18 12]); 
    subplot(1,2,1);  
       IMAGESC(x1+dx/2,x2,Rfl{1});
       title('\fontsize{13} R^{(1)}');
    subplot(1,2,2);  
       IMAGESC(x1+dx/2,x2, abs(dwE{1}));   
       title('\fontsize{13} abs(dwE^{(1)})');
       
    set(figure,'Units','centimeters','Position',[5 5 18 12]); 
    subplot(1,2,1);  
       IMAGESC(x1,x2+dx/2,Rfl{2});
       title('\fontsize{13} R^{(2)}');
    subplot(1,2,2);  
       IMAGESC(x1,x2+dx/2,abs(dwE{2}));   
       title('\fontsize{13} abs(dwE^{(2)})');
      
elseif nDIM == 3   
    N3cross = floor(input.N3/2+1);           % plot at x3 = 0 or x3 = dx/2
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1); 

    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,2,1);  
       IMAGESC(x1+dx/2,x2, Rfl{1}(:,:,N3cross));
       title('\fontsize{13} R^{(1)}');
    subplot(1,2,2);  
       IMAGESC(x1+dx/2,x2,abs(dwE{1}(:,:,N3cross)));
       title('\fontsize{13} abs(dwE^{(1)})');
       
    set(figure,'Units','centimeters','Position',[5 5 18 12]);         
    subplot(1,2,1);  
       IMAGESC(x1,x2+dx/2, Rfl{2}(:,:,N3cross));
       title('\fontsize{13} R^{(2)}');
    subplot(1,2,2);  
       IMAGESC(x1,x2+dx/2,abs(dwE{2}(:,:,N3cross)));
       title('\fontsize{13} abs(dwE^{(2)})');       
end; % if