function DisplayParameters
global Pnom P_it          
set(figure(9), 'Position', [0 0 800 500] );  N = length(P_it(:,1))-1;
subplot(1,4,1)
  plot(0:N,Pnom(1)*ones(N+1,1),':k','Linewidth',2);    hold on;
  plot(0:N,P_it(:,1),'-red','Linewidth',2);       
  title('P_1 = x^O_1','Fontsize',11); 
  xlabel('number of iterations','Fontsize',11); 
  axis('tight'); axis([ 0 N 0 25]); grid on; 
  legend('P nominal ', 'P iterated ','Location','SouthEast');
subplot(1,4,2)
  plot(0:N,Pnom(2)*ones(N+1,1),':k','Linewidth',2);    hold on;
  plot(0:N,P_it(:,2),'-red','Linewidth',2);      
  title('P_2 = x^O_2','Fontsize',11); 
  xlabel('number of iterations','Fontsize',11);
  axis('tight'); axis([ 0 N 0 25]);  grid on;
  legend('P nominal ', 'P iterated ','Location','SouthEast'); hold off;
subplot(1,4,3)
  plot(0:N,Pnom(3)*ones(N+1,1),':k','Linewidth',2);    hold on;
  plot(0:N,P_it(:,3),'-red','Linewidth',2);       
  title('P_3 = 1 - c_0^2/c_{sct}^2','Fontsize',11); 
  xlabel('number of iterations','Fontsize',11); 
  axis('tight'); axis([ 0 N -0.8 0.8]);  grid on;
  legend('P nominal ', 'P iterated ','Location','East');
subplot(1,4,4)
  plot(0:N,Pnom(4)*ones(N+1,1),':k','Linewidth',2);    hold on;
  plot(0:N,P_it(:,4),'-red','Linewidth',2);     
  title('P_4 = a','Fontsize',11); 
  xlabel('number of iterations','Fontsize',11); 
  axis('tight'); axis([ 0 N 20 80]);  grid on; 
  legend('P nominal ', 'P iterated ','Location','NorthEast')