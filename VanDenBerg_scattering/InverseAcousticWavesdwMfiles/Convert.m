function [] = Convert(Rfl,it,input)          
x1 = input.X1(:,1);   x2 = input.X2(1,:);    dx = input.dx;

MA = ones(1,2)/2;               % apply a two-points moving average to Rfl 
for j = 1:input.N2
   Rfl{1}(:,j)  = conv( Rfl{1}(1:input.N1,j),MA,'same'); 
   Rfl{2}(:,j)  = conv( Rfl{2}(1:input.N1,j),MA,'same'); 
end
for i = 1:input.N1
    Rfl{1}(i,:)  = conv( Rfl{1}(i,1:input.N2),MA,'same');  
    Rfl{2}(i,:)  = conv( Rfl{2}(i,1:input.N2),MA,'same'); 
end
                     
CHI_rho = RfltoCHI(Rfl,input);  % convert to mass-density volume contrast

set(figure,'Units','centimeters','Position',[0 2 17 9]); snapnow;

sp = subplot(1,3,1); 
     IMAGESC(x1+dx/2,x2,Rfl{1}); xlabel(''), ylabel('');
     title('\fontsize{12} R^{(1)}'); colorbar('hor');
     text(-70, -100,num2str(it),'FontSize',15);
     text(-100, -100,'it = ',    'FontSize',15);
     pos = get(sp,'Position');   pos = pos + [0.02 0 0 0]; 
     set(sp,'Position',pos);     
   
sp = subplot(1,3,2);
     IMAGESC(x1,x2+dx/2,Rfl{2});  xlabel(''), ylabel('');
     title('\fontsize{12} R^{(2)}');  colorbar('hor');  
     pos = get(sp,'Position');   pos = pos + [0 0 0 0]; 
     set(sp,'Position',pos);
      
sp = subplot(1,3,3);    
     IMAGESC(x1,x2,real(CHI_rho));      xlabel(''); ylabel('');
     title('\fontsize{13} 1-\rho_{sct}/ \rho_{0}');  colorbar('hor');
     pos = get(sp,'Position');   pos = pos - [0.02 0 0 0]; 
     set(sp,'Position',pos);