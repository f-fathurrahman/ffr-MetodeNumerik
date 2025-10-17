function [] = MovingAverage(Rfl,input)
x1  = input.X1(:,1);   x2 = input.X2(1,:); 

CHI = RfltoCHI(Rfl,input); 

set(figure,'Units','centimeters','Position',[0 -1 10 8]);  
    IMAGESC(x1,x2,CHI);  xlabel(''); ylabel('');  
set(figure,'Units','centimeters','Position',[0 0 10 8]);
    mesh(real(CHI),'edgecolor','k');
    view(37.5,45); axis xy; axis('off'); axis('tight');
    axis([1 input.N1 1 input.N2 [-1 1.01]*max(CHI(:))]);
        
  N = 2;  MA = ones(1,N)/N;   % two points average
  for k = 1:5
    for j = 1:input.N2
      CHI(:,j)  = conv(CHI(1:input.N1,j),MA,'same');  % moving average
    end
    for i = 1:input.N1
      CHI(i,:)  = conv(CHI(i,1:input.N2),MA,'same');  % moving average
    end
  end %k_loop
  
set(figure,'Units','centimeters','Position',[0 -1 10 8]);  
    IMAGESC(x1,x2,CHI);  xlabel(''); ylabel('');  
set(figure,'Units','centimeters','Position',[0 0 10 8]);
    mesh(real(CHI),'edgecolor','k');
    view(37.5,45); axis xy; axis('off'); axis('tight');
    axis([1 input.N1 1 input.N2 [-1 1.01]*max(CHI(:))]);
     