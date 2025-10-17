function v = gradj(p,j,input)
global nDIM;  dx = input.dx;  

if nDIM == 2;           
  N1 = input.N1;  N2 = input.N2;
    
  v = zeros(size(p));
  if j==1;  v(2:N1-1,:) = (p(3:N1,:) - p(1:N1-2,:)) / (2*dx); end;  % d1_p
  if j==2;  v(:,2:N2-1) = (p(:,3:N2) - p(:,1:N2-2)) / (2*dx); end;  % d2_p
             
elseif nDIM == 3;          
  N1 = input.N1;  N2 = input.N2;  N3 = input.N3;
    
  v = zeros(size(p)); 
  if j==1; v(2:N1-1,:,:)= (p(3:N1,:,:)-p(1:N1-2,:,:))/(2*dx); end;  % d1_p
  if j==2; v(:,2:N2-1,:)= (p(:,3:N2,:)-p(:,1:N2-2,:))/(2*dx); end;  % d2_p
  if j==3; v(:,:,2:N3-1)= (p(:,:,3:N3)-p(:,:,1:N3-2))/(2*dx); end;  % d3_p
end;