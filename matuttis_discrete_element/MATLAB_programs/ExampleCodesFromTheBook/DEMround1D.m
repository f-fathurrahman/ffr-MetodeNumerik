function [dydt]=DEMround1D(t,y);
% PURPOSE: Force computation function DEMround1D
%  (Program 7.5 from the book) 
% USAGE: To be called with the driver 
%  ch7p229roundparticlesimulation7_4.m
%  (Program 7.4 from the book)
% CAVEAT: If numerical dissipation should be introduced,
%  the non-smoothness should be dealt with as in
%  bouncing_ball_dissipation7_4.m
%  for point-particles
% REVISION HISTORY: 22-May-2014 H.-G. Matuttis

global m rad E lmax lmin g
n_part=length(m);
if length(y)~=2*length(m)
  error('length of y must be twice the length of m')
end
if length(rad)~=length(m)
  error('length of r must be the length of m')
end
a=zeros(1,n_part);
for i_part=1:n_part
  x1=y(2*i_part-1); % position of first particle
  rad1=rad(i_part);
  % Particle-Particle Interaction
  for j_part=i_part+1:n_part
    x2=y(2*j_part-1); % position of second particle
    rad2=rad(j_part);
    if (abs(x2-x1)<(rad(i_part)+rad(j_part))) % overlap
      forcemagnitude=E*abs(abs(x1-x2)-(rad1+rad2));
      forcedirection=sign(x1-x2);
      f=forcemagnitude*forcedirection;
      a(i_part)=a(i_part)+f;
      a(j_part)=a(j_part)-f; % use action=reaction
    end
  end
% Particle-wall Interaction
  if (x1-rad1)<lmin
    a(i_part)=a(i_part)-E*((x1-rad1)-lmin);
  end
  if (x1+rad1)>lmax
    a(i_part)=a(i_part)-E*((x1+rad1)-lmax);
  end
end
a=a+g;
dydt=zeros(2*n_part,1);
dydt(1:2:2*n_part-1)=y(2:2:2*n_part);
dydt(2:2:2*n_part)=a./m;
return