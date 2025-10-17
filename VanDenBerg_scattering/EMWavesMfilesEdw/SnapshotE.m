% This script is a continuation of TimedomainBiCGSTABFFTEdw 
clc; 

xS(1) = input.xS(1);  xS(2) = input.xS(2);
set(figure,'Units','centimeters','Position',[0 -1 18 32]);

for n =  1 : 12;
    
  subplot(4,3,n)
  n_t = (10 + n * 20);
  
% Compute the time domain power in dB
EdB = -20 * log10(Etime/max(max(Etime(:,:,n_t)))); 

% Make snapshot
  imagesc(input.X2(1,:),input.X1(:,1),EdB(:,:,n_t)); grid on; 
  t = n_t * input.dt *1e6;       % in microseconds
  title(['t = ',num2str(t),' {\mu}s']);  
  hold on;  colormap(hot);   axis('square');  caxis([0, 60]);
  plot(xS(2),xS(1),'ko','MarkerFaceColor',[0 0 0],'MarkerSize',3);
  phi =0:0.001:2*pi; 
  plot(input.a*cos(phi),input.a*sin(phi),'b','LineWidth',1.2);
  pause(0.1) 
  hold off;

end % if