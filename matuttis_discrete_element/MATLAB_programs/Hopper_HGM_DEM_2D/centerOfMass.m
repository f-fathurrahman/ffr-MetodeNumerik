% Calculate the center of mass of the particles.
% Walls are also treated as particles.
function [center_of_mass_x, center_of_mass_y, area] = centerOfMass(x, y)
% If if the corners x,y are not periodic, add a new one.
% modified HG Matuttis Thu Apr 11 20:20:59 JST 2013

% ???? modify if-condition
%
% Temporary center point for the particle.
temp_center_x=mean(x);
temp_center_y=mean(y);

% If the last corner is not the same corner as the first
% (as decided by the square of the
% relative distance to the average
% coordinate), the first corner is added to eliminate
% if-conditions due to the periodicity
tolerance=1e-10;
tolerance2=tolerance*tolerance;
veccenterx=temp_center_x-x;
veccentery=temp_center_y-y;
maxradius=max(veccenterx.^2+veccentery.^2);
distfirstend=(x(1)-x(end))^2+(y(1)-y(end))^2;

if (distfirstend>tolerance2*maxradius) 
  corner_num=length(x);
  x(corner_num+1)=x(1);
  y(corner_num+1)=y(1);
else
   corner_num=length(x)-1;
end
  

rx=x(1:corner_num)-temp_center_x;
ry=y(1:corner_num)-temp_center_y;
sx=x(2:corner_num+1)-x(1:corner_num);
sy=y(2:corner_num+1)-y(1:corner_num);

S=.5*(rx.*sy-ry.*sx);
cx=(x(1:corner_num)+x(2:corner_num+1)+temp_center_x)/3;
cy=(y(1:corner_num)+y(2:corner_num+1)+temp_center_y)/3;
 

center_of_mass_x = cx*S'/sum(S);
center_of_mass_y = cy*S'/sum(S);
area = sum(S);
return
