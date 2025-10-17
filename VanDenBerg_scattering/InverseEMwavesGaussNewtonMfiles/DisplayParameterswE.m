function DisplayParameterswE
global Pnom P_it       

set(figure(9), 'Position', [50 50 550 400] );  N = length(P_it(:,1))-1;

subplot(1,2,1)

  plot(0:N,Pnom(1)*ones(N+1,1),':k','Linewidth',2);    hold on;
  plot(0:N,P_it(:,1),'--b','Linewidth',2);       
  title('P_1 = eps_{sct}/eps_0','Fontsize',11); 
  xlabel('number of iterations','Fontsize',11); 
  axis('tight'); axis([ 0 N 1 1.7]);  grid on;

subplot(1,2,2)

  plot(0:N,Pnom(2)*ones(N+1,1),':k','Linewidth',2);    hold on;
  plot(0:N,P_it(:,2),'--b','Linewidth',2);     
  title('P_2 = a','Fontsize',11); 
  xlabel('number of iterations','Fontsize',11); 
  axis('tight'); axis([ 0 N 30 70]);  grid on; 
  legend('P nominal ', 'P wE ','Location','NorthEast')