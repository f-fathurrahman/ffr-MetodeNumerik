input = initAC;            x1 = input.X1(:,1);   x2 = input.X2(1,:);
close all

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
    title('\fontsize{12} it = 8');   axis('off'); axis(Axis);  view(View);
    pt=get(sp,'Position'); pt= pt+[0.065 0 0 0]; set(sp,'Position',pt);                 
sp=subplot(3,3,2);  mesh(CHI_16,'edgecolor', 'k');  
    title('\fontsize{12} it = 16');  axis('off'); axis(Axis);  view(View);
    pt=get(sp,'Position'); pt=pt+[0.00 0 0 0]; set(sp,'Position',pt);                 
sp=subplot(3,3,3);  mesh(CHI_32,'edgecolor', 'k');  
    title('\fontsize{12} it = 32');  axis('off'); axis(Axis);  view(View);
  pt=get(sp,'Position'); pt= pt+[-0.065 0 0 0]; set(sp,'Position',pt);

set(figure,'Units','centimeters','Position',[0 0 18 16]);
sp=subplot(3,3,1);  mesh(CHI_64,'edgecolor', 'k');  
    title('\fontsize{12} it = 64');  axis('off'); axis(Axis);  view(View);
    pt=get(sp,'Position'); pt=pt+[0.065 0 0 0]; set(sp,'Position',pt);                 
sp=subplot(3,3,2);  mesh(CHI_128,'edgecolor', 'k');  
    title('\fontsize{12} it = 128'); axis('off'); axis(Axis);  view(View);
    pt=get(sp,'Position'); pt=pt+[0.00 0 0 0]; set(sp,'Position',pt);                   
sp=subplot(3,3,3);  mesh(CHI_256,'edgecolor', 'k');  
    title('\fontsize{12} it = 256'); axis('off'); axis(Axis);  view(View);            
    pt=get(sp,'Position'); pt=pt+[-0.065 0 0 0]; set(sp,'Position',pt);                 