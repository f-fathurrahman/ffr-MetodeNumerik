close all
input = initAC; 

x1 = input.X1(:,1);   x2 = input.X2(1,:);   dx = input.dx;

load Contrast0   Rfl_0;    Convert(Rfl_0,0,input) ;  
load Contrast32  Rfl_32;   Convert(Rfl_32,32,input) ;    
load Contrast64  Rfl_64;   Convert(Rfl_64,64,input) ; 
load Contrast128 Rfl_128;  Convert(Rfl_128,128,input) ; 
load Contrast256 Rfl_256;  Convert(Rfl_256,256,input) ; 

% Plot exact case
set(figure,'Units','centimeters','Position',[0 2 17 9]); snapnow;

sp = subplot(1,3,1); 
     IMAGESC(x1+dx/2,x2,real(input.Rfl{1}));  xlabel(''); ylabel('');
     title('\fontsize{12} R^{(1)}');                     colorbar('hor');
     text(-100, -100,'Exact ',    'FontSize',15);
     pos = get(sp,'Position');   pos = pos + [0.02 0 0 0]; 
     set(sp,'Position',pos);     
     
sp = subplot(1,3,2);
     IMAGESC(x1,x2+dx/2,real(input.Rfl{2}));  xlabel(''); ylabel('');
     title('\fontsize{12} R^{(2)}');                     colorbar('hor');
     pos = get(sp,'Position');   pos = pos + [0 0 0 0];
     set(sp,'Position',pos);
   
sp = subplot(1,3,3);    
     IMAGESC(x1,x2,real(input.CHI_rho));      xlabel(''); ylabel('');
     title('\fontsize{13} 1-\rho_{sct}/ \rho_{0}');      colorbar('hor');
     pos = get(sp,'Position');   pos = pos - [0.02 0 0 0]; 
     set(sp,'Position',pos);