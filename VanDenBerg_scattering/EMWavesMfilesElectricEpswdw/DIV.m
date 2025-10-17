function u = DIV(v,input)
global nDIM;  dx = input.dx;  

if nDIM == 1;               
  N1 = input.N1;
    
  u = zeros(size(v{1})); 
  u(2:N1-1) = (v{1}(3:N1) - v{1}(1:N1-2))/(2*dx);                  % d1v1
   
elseif nDIM == 2;         
  N1 = input.N1;  N2 = input.N2;
    
  v{1}(2:N1-1,:) =  v{1}(3:N1,:) - v{1}(1:N1-2,:);                 % d1v1
  v{2}(:,2:N2-1) =  v{2}(:,3:N2) - v{2}(:,1:N2-2);                 % d2v2
               u = (v{1} + v{2}) / (2*dx);

elseif nDIM == 3;           
  N1 = input.N1;  N2 = input.N2;  N3 = input.N3;
    
  v{1}(2:N1-1,:,:) = v{1}(3:N1,:,:) - v{1}(1:N1-2,:,:);            % d1v1
  v{2}(:,2:N2-1,:) = v{2}(:,3:N2,:) - v{2}(:,1:N2-2,:);            % d2v2
  v{3}(:,:,2:N3-1) = v{3}(:,:,3:N3) - v{3}(:,:,1:N3-2);            % d3v3
                   u = (v{1} + v{2} + v{3}) / (2*dx);
end;
