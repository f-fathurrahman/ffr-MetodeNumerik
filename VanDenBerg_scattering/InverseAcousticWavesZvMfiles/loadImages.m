input = initAC;    x1 = input.X1(:,1);   x2 = input.X2(1,:);
close all

set(figure,'Units','centimeters','Position',[0 -1 20 30]);  
sp = subplot(3,3,1);   load Contrast0 CHI_0;    IMAGESC(x1,x2,CHI_0);   
     title('\fontsize{13} it = 0');  xlabel(''); ylabel(''); 
     pos=get(sp,'Position'); pos=pos+[0.05 0 0 0]; set(sp,'Position',pos);
sp = subplot(3,3,2);   load Contrast1   CHI_1;  IMAGESC(x1,x2,CHI_1);   
     title('\fontsize{13} it = 1');  xlabel(''); ylabel(''); 
     pos=get(sp,'Position'); pos=pos+[0 0 0 0];    set(sp,'Position',pos);   
sp = subplot(3,3,3);   load Contrast4   CHI_4;  IMAGESC(x1,x2,CHI_4);  
     title('\fontsize{13} it = 4');  xlabel(''); ylabel(''); 
     pos=get(sp,'Position');pos=pos+[-0.05 0 0 0]; set(sp,'Position',pos); 
  
set(figure,'Units','centimeters','Position',[0 -1 20 30]);   
sp = subplot(3,3,1);   load Contrast8   CHI_8;  IMAGESC(x1,x2,CHI_8);  
     title('\fontsize{13} it = 8');  xlabel(''); ylabel(''); 
     pos=get(sp,'Position'); pos=pos+[0.05 0 0 0]; set(sp,'Position',pos);   
sp = subplot(3,3,2);   load Contrast16  CHI_16; IMAGESC(x1,x2,CHI_16);
     title('\fontsize{13} it = 16');  xlabel(''); ylabel('');  
     pos=get(sp,'Position'); pos=pos+[0 0 0 0];    set(sp,'Position',pos);  
sp = subplot(3,3,3);   load Contrast32  CHI_32; IMAGESC(x1,x2,CHI_32);  
     title('\fontsize{13} it = 32');  xlabel(''); ylabel(''); 
     pos=get(sp,'Position');pos=pos+[-0.05 0 0 0]; set(sp,'Position',pos);

set(figure,'Units','centimeters','Position',[0 -1 20 30]);    
sp = subplot(3,3,1);  load Contrast64  CHI_64;  IMAGESC(x1,x2,CHI_64);  
   title('\fontsize{13} it = 64');  xlabel(''); ylabel('');  
     pos=get(sp,'Position'); pos=pos+[0.05 0 0 0]; set(sp,'Position',pos); 
sp = subplot(3,3,2);  load Contrast128 CHI_128; IMAGESC(x1,x2,CHI_128); 
   title('\fontsize{13} it = 128'); xlabel(''); ylabel('');
     pos=get(sp,'Position'); pos=pos+[0 0 0 0];    set(sp,'Position',pos); 
sp = subplot(3,3,3);  load Contrast256 CHI_256; IMAGESC(x1,x2,CHI_256); 
   title('\fontsize{13} it = 256'); xlabel(''); ylabel('');  
     pos=get(sp,'Position');pos=pos+[-0.05 0 0 0]; set(sp,'Position',pos);