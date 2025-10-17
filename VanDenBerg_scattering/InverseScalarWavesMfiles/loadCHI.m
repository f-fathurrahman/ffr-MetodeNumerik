input = init;            x1 = input.X1(:,1);   x2 = input.X2(1,:);

close all

set(figure,'Units','centimeters','Position',[0 -1 20 30]);   
sp = subplot(3,3,1);            load Contrast0 CHI_0;    
   IMAGESC(x1,x2,CHI_0);   title('\fontsize{13} it = 0');  xlabel(''); ylabel(''); 
   pos=get(sp,'Position'); pos=pos + [0.05 0 0 0]; set(sp,'Position',pos);
sp =subplot(3,3,2);            load Contrast1   CHI_1;    
   IMAGESC(x1,x2,CHI_1);   title('\fontsize{13} it = 1');  xlabel(''); ylabel(''); 
   pos=get(sp,'Position'); pos=pos +[0 0 0 0]; set(sp,'Position',pos);   
sp = subplot(3,3,3);            load Contrast4   CHI_4;
   IMAGESC(x1,x2,CHI_4);   title('\fontsize{13} it = 4');  xlabel(''); ylabel(''); 
       pos=get(sp,'Position'); pos=pos +[-0.05 0 0 0]; set(sp,'Position',pos); 
  
set(figure,'Units','centimeters','Position',[0 -1 20 30]);   
sp =subplot(3,3,1);            load Contrast8   CHI_8;    
   IMAGESC(x1,x2,CHI_8);   title('\fontsize{13} it = 8');  xlabel(''); ylabel(''); 
     pos=get(sp,'Position'); pos=pos +[0.05 0 0 0]; set(sp,'Position',pos);   
sp = subplot(3,3,2);            load Contrast16  CHI_16;   
   IMAGESC(x1,x2,CHI_16);  title('\fontsize{13} it = 16');  xlabel(''); ylabel('');  
    pos=get(sp,'Position'); pos=pos +[0 0 0 0]; set(sp,'Position',pos);  
sp = subplot(3,3,3);            load Contrast32  CHI_32;   
   IMAGESC(x1,x2,CHI_32);  title('\fontsize{13} it = 32');  xlabel(''); ylabel(''); 
     pos=get(sp,'Position'); pos=pos +[-0.05 0 0 0]; set(sp,'Position',pos);

set(figure,'Units','centimeters','Position',[0 -1 20 30]);    
sp = subplot(3,3,1);            load Contrast64  CHI_64;   
   IMAGESC(x1,x2,CHI_64);  title('\fontsize{13} it = 64');  xlabel(''); ylabel('');  
     pos=get(sp,'Position'); pos=pos +[0.05 0 0 0]; set(sp,'Position',pos); 
sp = subplot(3,3,2);            load Contrast128 CHI_128;  
   IMAGESC(x1,x2,CHI_128); title('\fontsize{13} it = 128'); xlabel(''); ylabel('');
     pos=get(sp,'Position'); pos=pos +[0 0 0 0]; set(sp,'Position',pos); 
sp = subplot(3,3,3);            load Contrast256 CHI_256;  
   IMAGESC(x1,x2,CHI_256); title('\fontsize{13} it = 256'); xlabel(''); ylabel('');  
     pos=get(sp,'Position'); pos=pos +[-0.05 0 0 0]; set(sp,'Position',pos);
  


Axis = [1 input.N1 1 input.N2 [-1 1.01]*max(CHI_256(:))]; View=[37.5,45]; 

set(figure,'Units','centimeters','Position',[0 0 19 16]);
sp = subplot(3,3,1);  mesh(CHI_0,'edgecolor', 'k');   
    title('\fontsize{12} it = 0');   axis('off'); axis(Axis); view(View);
    pt=get(sp,'Position'); pt=pt +[0.065 0 0 0]; set(sp,'Position',pt);                 
sp = subplot(3,3,2);  mesh(CHI_1,'edgecolor', 'k');  
    title('\fontsize{12} it = 1');   axis('off'); axis(Axis); view(View);
    pt=get(sp,'Position'); pt=pt +[0.00 0 0 0]; set(sp,'Position',pt);                 
sp = subplot(3,3,3);  mesh(CHI_4,'edgecolor', 'k');   
    title('\fontsize{12} it = 4');   axis('off'); axis(Axis); view(View);
    pt=get(sp,'Position'); pt=pt+[-0.065 0 0 0]; set(sp,'Position',pt); 

set(figure,'Units','centimeters','Position',[0 0 18 16]); 
sp=subplot(3,3,1);  mesh(CHI_8,'edgecolor', 'k');  
    title('\fontsize{12} it = 8');   axis('off'); axis(Axis);view(View);
    pt=get(sp,'Position'); pt= pt+[0.065 0 0 0]; set(sp,'Position',pt);                 
sp=subplot(3,3,2);  mesh(CHI_16,'edgecolor', 'k');  
    title('\fontsize{12} it = 16');  axis('off'); axis(Axis);  view(View);
    pt=get(sp,'Position'); pt=pt+[0.00 0 0 0]; set(sp,'Position',pt);                 
sp=subplot(3,3,3);  mesh(CHI_32,'edgecolor', 'k');  
    title('\fontsize{12} it = 32');  axis('off'); axis(Axis); view(View);
  pt=get(sp,'Position'); pt= pt+[-0.065 0 0 0]; set(sp,'Position',pt);

set(figure,'Units','centimeters','Position',[0 0 18 16]);
sp=subplot(3,3,1);  mesh(CHI_64,'edgecolor', 'k');  
    title('\fontsize{12} it = 64');  axis('off'); axis(Axis); view(View);
    pt=get(sp,'Position'); pt=pt+[0.065 0 0 0]; set(sp,'Position',pt);                 
sp=subplot(3,3,2);  mesh(CHI_128,'edgecolor', 'k');  
    title('\fontsize{12} it = 128'); axis('off'); axis(Axis); view(View);
    pt=get(sp,'Position'); pt=pt+[0.00 0 0 0]; set(sp,'Position',pt);                   
sp=subplot(3,3,3);  mesh(CHI_256,'edgecolor', 'k');  
    title('\fontsize{12} it = 256'); axis('off'); axis(Axis); view(View);            
    pt=get(sp,'Position'); pt=pt+[-0.065 0 0 0]; set(sp,'Position',pt);                 
                 