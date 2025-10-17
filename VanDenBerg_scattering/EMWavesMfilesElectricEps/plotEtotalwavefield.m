function  plotEtotalwavefield(E_inc,E_sct,input)
global nDIM;   
a = input.a;    

    E = cell(1,nDIM);
for n = 1:nDIM 
    E{n} = E_inc{n} + E_sct{n};
end

if nDIM == 2   % Plot wave fields in two-dimensional space ----------------

   set(figure,'Units','centimeters','Position',[5 5 18 12]);
   x1 = input.X1(:,1);   x2 = input.X2(1,:); 
   N1 = input.N1;        N2 = input.N2;
   subplot(1,2,1); 
      IMAGESC(x1,x2, abs(E{1})); 
      title(['\fontsize{13} 2D Electric field E_1 '])         
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');
   subplot(1,2,2);
      IMAGESC(x1,x2, abs(E{2})); 
      title(['\fontsize{13} 2D Electric field E_2 '])            
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');

elseif nDIM == 3   % Plot wave fields at x3 = 0 or x3 = dx/2 --------------
    
   set(figure,'Units','centimeters','Position',[5 5 18 12]);                                    
    N1 = input.N1;         N2 = input.N2;    N3cross = floor(input.N3/2+1);
    x1 = input.X1(:,1,1);  x2 = input.X2(1,:,1);   
   subplot(1,3,1)
      IMAGESC(x1,x2,abs(E{1}(:,:,N3cross)));
      title(['\fontsize{13} 3D Electric field E_1 '])     
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w');    
   subplot(1,3,2)    
      IMAGESC(x1(2:N1-1),x2(2:N2-1),abs(E{2}(2:N1-1,2:N2-1,N3cross)));
      title(['\fontsize{13} 3D Electric field E_2 ']); 
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w'); 
   subplot(1,3,3)    
      IMAGESC(x1(2:N1-1),x2(2:N2-1),abs(E{3}(2:N1-1,2:N2-1,N3cross)));
      title(['\fontsize{13} 3D Electric field E_3 ']); 
      hold on; phi = 0:.01:2*pi; plot(a*cos(phi),a*sin(phi),'w'); 
end % if
