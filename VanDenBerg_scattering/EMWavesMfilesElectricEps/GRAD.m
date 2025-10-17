function v = GRAD(p,input)
global nDIM;  dx = input.dx;  

if nDIM == 1;              N1 = input.N1;
    
  v{1} = zeros(size(p));
  v{1}(2:N1-1) = (p(3:N1) - p(1:N1-2)) / (2*dx);                    % d1_p
  
elseif nDIM == 2;          N1 = input.N1;  N2 = input.N2;
    
  v{1} = zeros(size(p)); v{2} = v{1};
  v{1}(2:N1-1,:) = (p(3:N1,:) - p(1:N1-2,:)) / (2*dx);              % d1_p
  v{2}(:,2:N2-1) = (p(:,3:N2) - p(:,1:N2-2)) / (2*dx);              % d2_p
             
elseif nDIM == 3;          N1 = input.N1;  N2 = input.N2;  N3 = input.N3;
    
  v{1} = zeros(size(p)); v{2} = v{1};   v{3} = v{1};
  v{1}(2:N1-1,:,:) = (p(3:N1,:,:) - p(1:N1-2,:,:)) / (2*dx);        % d1_p
  v{2}(:,2:N2-1,:) = (p(:,3:N2,:) - p(:,1:N2-2,:)) / (2*dx);        % d2_p
  v{3}(:,:,2:N3-1) = (p(:,:,3:N3) - p(:,:,1:N3-2)) / (2*dx);        % d3_p
end