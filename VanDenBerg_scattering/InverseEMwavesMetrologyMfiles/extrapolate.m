function [u] = extrapolate(v,grid_out,grid_in,input)
global nDIM;                                  u = zeros(size(v));
 
if nDIM == 1;        N1 = input.N1; 
    
  if grid_out == 1;  %--------------------------------forward extrapolation
     u(2:N1-2) = v(2:N1-2) + (v(3:N1-1)-v(1:N1-3))/4;
  elseif grid_in ==1 %-------------------------------backward extrapolation
     u(2:N1-2) = v(2:N1-2) - (v(3:N1-1)-v(1:N1-3))/4; 
  end;
  
elseif nDIM == 2;    N1 = input.N1;  N2 = input.N2; 
    
  if grid_out == 1;  %--------------------------------forward extrapolation
     u(2:N1-2,:) = v(2:N1-2,:) + (v(3:N1-1,:)-v(1:N1-3,:))/4;    
  elseif grid_out == 2;  
    u(:,2:N2-2) = v(:,2:N2-2) + (v(:,3:N2-1)-v(:,1:N2-3))/4;               
  end;          
  if grid_in == 1;  %--------------------------------backward extrapolation
      u(2:N1-2,:) = v(2:N1-2,:) - (v(3:N1-1,:)-v(1:N1-3,:))/4; 
  elseif grid_in == 2; 
     u(:,2:N2-2) = v(:,2:N2-2) - (v(:,3:N2-1)-v(:,1:N2-3))/4;
  end;          
 
elseif nDIM == 3;   N1 = input.N1;  N2 = input.N2;  N3 = input.N3;  
    
  if grid_out == 1; %---------------------------------forward extrapolation
     u(2:N1-1,:,:) = v(2:N1-1,:,:) + (v(3:N1,:,:)-v(1:N1-2,:,:))/4;
  elseif grid_out == 2;                         
     u(:,2:N2-1,:) = v(:,2:N2-1,:) + (v(:,3:N2,:)-v(:,1:N2-2,:))/4;
  elseif grid_out == 3;                         
     u(:,:,2:N3-1) = v(:,:,2:N3-1) + (v(:,:,3:N3)-v(:,:,1:N3-2))/4;
  end;          
  if grid_in == 1;  %--------------------------------backward extrapolation
    u(3:N1-2,:,:) = v(3:N1-2,:,:) - (v(4:N1-1,:,:)-v(2:N1-3,:,:))/4;
  elseif grid_in == 2; 
    u(:,3:N2-2,:) = v(:,3:N2-2,:) - (v(:,4:N2-1,:)-v(:,2:N2-3,:))/4;         
  elseif grid_in == 3; 
    u(:,:,3:N3-2) = v(:,:,3:N3-2) - (v(:,:,4:N3-1)-v(:,:,2:N3-3))/4;         
  end;           
end;   