function [dydt]=DEMround2Dnorot(t,y);
% USAGE: to be called from
%  ch7p232roundparticlesimulation7_6
% CAVEAT: contains no neighborhood-algorithm
%  and computes the interaction "brute force"
%  with two nested loops which may take
%  1/2*n_part^2
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

global m rad E lmax lmin lmaxx lminx lmaxy lminy g

n_part=length(m);
if length(y)~=4*length(m)
  error('length of y must be four times the lenght of m')
end    
if length(rad)~=length(m)
  error('length of rad must be the lenght of m')
end    

a=zeros(2,n_part);
for i_part=1:n_part
  r1=[y(4*i_part-3)
      y(4*i_part-1)]; %  position of first particle
  rad1=rad(i_part);
% Particle-Particle Interaction
  for j_part=i_part+1:n_part 
    r2=[y(4*j_part-3)
        y(4*j_part-1)]; % position of second particle
    rad2=rad(j_part);
    if (norm(r1-r2)<(rad(i_part)+rad(j_part)))
      forcemagnitude=E*abs(norm(r1-r2)-(rad1+rad2));
      forcedirection=(r1-r2)/norm(r1-r2);  
      f=forcemagnitude*forcedirection;
      a(:,i_part)=a(:,i_part)+f;
      a(:,j_part)=a(:,j_part)-f;
    end    
  end
% Particle-wall Interaction  
  if (r1(1)-rad1)<lminx
    a(1,i_part)=a(1,i_part)-E*((r1(1)-rad1)-lminx);  
  end    
  if (r1(1)+rad1)>lmaxx
    a(1,i_part)=a(1,i_part)-E*((r1(1)+rad1)-lmaxx);  
  end    
  if (r1(2)-rad1)<lminy
    a(2,i_part)=a(2,i_part)-E*((r1(2)-rad1)-lminy);  
  end    
  if (r1(2)+rad1)>lmaxy
    a(2,i_part)=a(2,i_part)-E*((r1(2)+rad1)-lmaxy);  
  end    
end
a(2,:)=a(2,:)+g;
dydt=zeros(4*n_part,1);
dydt(1:4:4*n_part-3)=y(2:4:4*n_part-2); 
dydt(2:4:4*n_part-2)=a(1,:)./m;
dydt(3:4:4*n_part-1)=y(4:4:4*n_part); 
dydt(4:4:4*n_part)=a(2,:)./m;

return
