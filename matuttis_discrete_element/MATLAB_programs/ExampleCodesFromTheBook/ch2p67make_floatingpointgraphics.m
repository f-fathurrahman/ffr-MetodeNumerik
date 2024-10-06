clear all
format compact
clf
% PURPOSE: plot the distribution of floating-point
%  numbers with base=2, a mantissa of 3 and an 
%  exponent from -1 to 3, and below the distribution
%  of integer numbers in the same interval
% USAGE: Was used to produce Fig. 2.1
% ALGORITHM:  see p. 66
% REVISION HISTORY: 21-May-2014 H.-G. Matuttis


j=0
for i1=0:1 % first digit of the mantissa
  for i2=0:1 % second  digit of the mantissa
    for i3=0:1 % third digit of the mantissa
% to extend to larger mantissas, include
% additional loops
      for r=-1:3 % exponent
        j=j+1;
        mantissa=(i1/2+i2/4+i3/8);
        exponent=2^r;
        x1(j)=mantissa*exponent;
      end
    end
  end
end
x1=[-x1 x1];

% plot the floating-point values as vertical lines
plot([x1 ;x1],[0*x1 ;.5+0*x1],'k'  )
hold on
% plot the integers as vertical lines below
x2=[-7:7];
plot([x2 ;x2],[1.25+0*x2 ;1.75+0*x2],'k'  )
axis off
% annotate the integers
for i=min(x2):max(x2)
  text(i-.15+.051*sign(i+.5),2.1,num2str(i))  
end    

axis image

return
