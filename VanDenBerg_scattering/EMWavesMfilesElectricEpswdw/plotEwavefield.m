function  plotEwavefield(Title,E_inc,E,input)
global nDIM;   
a = input.a;        E_sct = E - E_inc;

if nDIM == 2;  % Plot wave fields in two-dimensional space ------------

   set(figure,'Units','centimeters','Position',[5 5 18 12]);
   x1 = input.X1(:,1);   x2 = input.X2(1,:); 
   N1 = input.N1;        N2 = input.N2;
   subplot(1,2,1); 
      IMAGESC(x1(2:N1-1),x2(2:N2-1), abs(E_sct(2:N1-1,2:N2-1))); 
      title(['\fontsize{13} 2D Scattered field: ',Title])         
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');
   subplot(1,2,2);
      IMAGESC(x1(2:N1-1),x2(2:N2-1), abs(E(2:N1-1,2:N2-1))); 
         title(['\fontsize{13} 2D Total field: ',Title])            
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');

elseif nDIM == 3;  % Plot wave fields at x3 = 0 or x3 = dx/2 --------------
    
   set(figure,'Units','centimeters','Position',[5 5 18 12]);                                    
    N1 = input.N1;         N2 = input.N2;    N3cross = floor(input.N3/2+1);
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1);   
   subplot(1,2,1)
      IMAGESC(x1(2:N1-1),x2(2:N2-1),abs(E_sct(2:N1-1,2:N2-1,N3cross)));
      title(['\fontsize{13} 3D Scattered field: ',Title])     
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');    
   subplot(1,2,2)    
      IMAGESC(x1(2:N1-1),x2(2:N2-1),abs(E(2:N1-1,2:N2-1,N3cross)));
      title(['\fontsize{13} 3D Total field: ',Title]); 
      caxis([ 0 max(abs(E_sct(:)))]);
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w'); 
      
end % if
