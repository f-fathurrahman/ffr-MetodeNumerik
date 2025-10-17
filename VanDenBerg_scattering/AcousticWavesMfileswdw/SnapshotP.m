% This script is a continuation of TimedomainBiCGSTABFFTwdw 
clc; 
xS(1) = input.xS(1);  xS(2) = input.xS(2);

ptimeabs = abs(ptime);
set(figure,'Units','centimeters','Position',[0 -1 18 32]);
for n =  1 : 12;
  subplot(4,3,n)
 n_t = (20 + n * 7);
  
% Compute the time domain power in dB
  udB = -20 * log10(ptimeabs/max(max(ptimeabs(:,:,n_t))));    
% Make snapshot
  imagesc(input.X2(1,:),input.X1(:,1),udB(:,:,n_t)); grid on; 
  title(['t = ',num2str(n_t * input.dt),' s']);  
  hold on;  colormap(hot);  axis('square');  caxis([0, 40]);
  plot(xS(2),xS(1),'ko','MarkerFaceColor',[0 0 0],'MarkerSize',3);
  phi = 0:0.001:2*pi; 
  plot(input.a*cos(phi),input.a*sin(phi),'b','LineWidth',1.2);
  pause(0.1)
  hold off; 

end % if