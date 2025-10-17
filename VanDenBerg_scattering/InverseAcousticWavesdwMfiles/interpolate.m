function [v] = interpolate(v,grid_out,grid_in,input)
global nDIM;  
   
if nDIM == 2;            N1 = input.N1;  N2 = input.N2; 
   
   if grid_in     == 1;  v(2:N1-1,:) = (v(1:N1-2,:) + v(2:N1-1,:))/2;    
   elseif grid_in == 2;  v(:,2:N2-1) = (v(:,1:N2-2) + v(:,2:N2-1))/2;                       
   end
             
   if grid_out     == 1; v(1:N1-1,:) = (v(1:N1-1,:) + v(2:N1,:))/2;  
   elseif grid_out == 2; v(:,1:N2-1) = (v(:,1:N2-1) + v(:,2:N2))/2;          
   end     
   
elseif nDIM == 3;        N1 = input.N1;  N2 = input.N2;  N3 = input.N3; 
  
   if grid_in     == 1;  v(2:N1-1,:,:) = (v(1:N1-2,:,:) + v(2:N1-1,:,:))/2;                               
   elseif grid_in == 2;  v(:,2:N2-1,:) = (v(:,1:N2-2,:) + v(:,2:N2-1,:))/2;                                                 
   elseif grid_in == 3;  v(:,:,2:N3-1) = (v(:,:,1:N3-2) + v(:,:,2:N3-1))/2;           
   end
    
   if grid_out     == 1; v(1:N1-1,:,:) = (v(1:N1-1,:,:) + v(2:N1,:,:))/2;      
   elseif grid_out == 2; v(:,1:N2-1,:) = (v(:,1:N2-1,:) + v(:,2:N2,:))/2;                       
   elseif grid_out == 3; v(:,:,1:N3-1) = (v(:,:,1:N3-1) + v(:,:,2:N3))/2;
   end  
end   

