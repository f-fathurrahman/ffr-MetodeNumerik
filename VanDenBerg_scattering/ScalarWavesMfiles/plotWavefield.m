function  plotWavefield(u_inc,u,input)
global nDIM;  
a = input.a;       u_sct = u - u_inc;

if nDIM == 1       % Plot wave fields in one-dimensional space ------------

   set(figure,'Units','centimeters','Position',[5 5 18 12]);
   x1 = input.X1;  N1 = input.N1;
   x2 = -30 : 30;  
   subplot(1,2,1)   
      IMAGESC(x1(2:N1-1),x2, abs(u_sct(2:N1-1)) );  
      title('\fontsize{13} 1D: abs(u^{sct})');  
      hold on;     plot(3*x2,0*x2+a,'--w','LineWidth',2); 
                   plot(3*x2,0*x2-a,'--w','LineWidth',2);
  subplot(1,2,2)    
      IMAGESC(x1(2:N1-1),x2, abs(u(2:N1-1)) );  
      title('\fontsize{13} 1D: abs(u)'); caxis([ 0 max(abs(u_sct))]);  
      hold on;     plot(3*x2,0*x2+a,'--w','LineWidth',2); 
                   plot(3*x2,0*x2-a,'--w','LineWidth',2); 

elseif nDIM == 2  % Plot wave fields in two-dimensional space -------------

    set(figure,'Units','centimeters','Position',[5 5 18 12]);
   x1 = input.X1(:,1);   x2 = input.X2(1,:); 
   N1 = input.N1;        N2 = input.N2;
   subplot(1,2,1); 
      IMAGESC(x1(2:N1-1),x2(2:N2-1), abs(u_sct(2:N1-1,2:N2-1))); 
      title('\fontsize{12} 2D: abs(u^{sct})');         
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');
   subplot(1,2,2);
      IMAGESC(x1(2:N1-1),x2(2:N2-1), abs(u(2:N1-1,2:N2-1))); 
      title('\fontsize{13} 2D: abs(u)'); caxis([ 0 max(abs(u_sct(:)))]);         
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');

elseif nDIM == 3  % Plot wave fields at x3 = 0 or x3 = dx/2 ---------------

   set(figure,'Units','centimeters','Position',[5 5 18 12]);                                    
    N1 = input.N1;         N2 = input.N2;         N3 = input.N3; 
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1);  N3cross = floor(N3/2+1); 
   subplot(1,2,1)
      IMAGESC(x1(2:N1-1),x2(2:N2-1),abs(u_sct(2:N1-1,2:N2-1,N3cross)));
      title('\fontsize{13} 3D:abs(u^{sct})'); 
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');    
   subplot(1,2,2)    
      IMAGESC(x1(2:N1-1),x2(2:N2-1),abs(u(2:N1-1,2:N2-1,N3cross)));
      title('\fontsize{13} 3D: abs(u)'); caxis([ 0 max(abs(u_sct(:)))]);
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');
 
end % if