function plotFieldError(Title,Field_exact,Field,x1,x2,input)  
global nDIM;  a = input.a;   PHI = 0:.01:2*pi;


if nDIM == 3;  N3cross =floor(input.N3/2+1); 
  Field_exact = Field_exact(:,:,N3cross);
  Field       = Field(:,:,N3cross);
end

set(figure,'Units','centimeters','Position',[5 5 18 10]);            
subplot(1,3,1)     
IMAGESC(x1,x2,abs(Field_exact)) 
   title(['\fontsize{10}', Title,'_{exact}'])
   hold on;  plot(a*cos(PHI),a*sin(PHI),'w');
subplot(1,3,2)   
IMAGESC(x1,x2,abs(Field))
   title(['\fontsize{10}', Title])
   hold on;  plot(a*cos(PHI),a*sin(PHI),'w'); 
subplot(1,3,3)
IMAGESC(x1,x2,abs(Field - Field_exact))
   title('\fontsize{10}  abs(Error)');   
   hold on; plot(a*cos(PHI),a*sin(PHI),'w');
     
