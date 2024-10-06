function [m,mi,f,x,y]= mass_momentinertia(...
    nsides,polx,poly,rho)
% [m,mi,mr,mr2,f,x,y]=mass_momentinertia(n,nsides,polx,poly, rho)
% PURPOSE: Calculates the masses, moments of inertia 
%  and centers of mass for a vector polygon-data 
%  from the area and the density.
% INPUTS:
%  nsides - Vector of number of corners for each particles.
%  polx - X-coordinates in each particles.
%  poly - Y-coordinates in each particles.
%  rho - Density of the particle.
% ALGORITHM: See below
% OUTPUTS:
%  m - Mass
%  mi - Moment of inertia with respect to the center of mass
%  f - Angle ("phi")
%  x - Centers of mass, x-coordinate
%  y - Centers of mass, y-coordinate
%
% RECORD OF REVISION:
%  12-May-2014   HG Simplified interfaces 
%  2013-Jun-11   S. H. Ng Translated from Fortran
%

n=length(nsides);
m =zeros(1,n);
mi=zeros(1,n);
f =zeros(1,n); 
x =zeros(1,n);
y =zeros(1,n);

onethird=1.d0/3.d0;
one_twelfth=1.d0/12.d0;

% Add periodicity to the x/y-coordinates of the corners
% if necessary 
tolerance=1e-10;
tolerance2=tolerance*tolerance;
for j=1:n
% Temporary center point for the particle.
  temp_center_x=mean(polx(1:nsides(j),j));
  temp_center_y=mean(poly(1:nsides(j),j));
% If the last corner is not the same corner as the first
% (as decided by the square of the
% relative distance to the average
% coordinate), the first corner is added to eliminate
% if-conditions due to the periodicity
  veccenterx=temp_center_x-polx(1:nsides(j),j);
  veccentery=temp_center_y-poly(1:nsides(j),j);
  maxradius=max(veccenterx.^2+veccentery.^2);
  distfirstend=(polx(1,j)-polx(nsides(j),j))^2+...
               (poly(1,j)-poly(nsides(j),j))^2;

  if (distfirstend>tolerance2*maxradius) 
    polx(nsides(j)+1,j)=polx(1,j);
    poly(nsides(j)+1,j)=poly(1,j);
  end
end

for j = 1:n
  area=0.0d0;
  sx  =0.0d0;
  sy  =0.0d0;
  for i= 2:nsides(j)-1
    x1=(polx(i,j)  -polx(1,j));
    x2=(polx(i+1,j)-polx(1,j));
    y1=(poly(i,j)  -poly(1,j));
    y2=(poly(i+1,j)-poly(1,j));
    ai=x1*y2 - x2*y1;
    sx=sx + ai*(x1+x2);
    sy=sy + ai*(y1+y2);
    area = area + ai;
  end
  sx = (1d0/3d0)*sx/area + polx(1,j);
  sy = (1d0/3d0)*sy/area + poly(1,j);
  area = 0.5d0*area;
   
  x(j) = sx;
  y(j) = sy;
   
  m(j) = area*rho;
end

% Formula for the momentum of inertia of a triangle:
% I=m/2*(s*s+l*l/12)
% where m is the mass, 
% s the length of the line which cuts the third
% edge in two equal parts and 
% l the length of this third edge.
% L2=l*l=(a*a+b*b-2ab cos fi)
% S2=s*s=(a*a+b*b+2ab cos fi) / 4
for j=1:n   
  mi(j)=0.d0;
  for i=2:nsides(j)+1 
    rayx =polx(i,j)-x(j);
    rayy =poly(i,j)-y(j);
    sidex=polx(i,j)-polx(i-1,j);
    sidey=poly(i,j)-poly(i-1,j);
    area_triangle=.5*(rayx*sidey-rayy*sidex);
    center_of_mass_triangle_x=onethird*((polx(i,j)-x(j))+(polx(i-1,j)-x(j)));
    center_of_mass_triangle_y=onethird*((poly(i,j)-y(j))+(poly(i-1,j)-y(j)));
% Composition of the momenta of inertia
% via Steiners Theorem
    center_of_mass_distance=...
    sqrt(center_of_mass_triangle_x^2+center_of_mass_triangle_y^2);
    delta_I_steiner=area_triangle*center_of_mass_distance^2;
    delta_I_own=area_triangle*one_twelfth*...
      (center_of_mass_distance^2+...
      (center_of_mass_triangle_x-(polx(i,j)-x(j)) )^2+...
      (center_of_mass_triangle_y-(poly(i,j)-y(j)) )^2+...
      (center_of_mass_triangle_x-(polx(i-1,j)-x(j)) )^2+...
      (center_of_mass_triangle_y-(poly(i-1,j)-y(j)) )^2 );
    mi(j)=mi(j)+delta_I_steiner+delta_I_own;      
  end
   
  mi(j)=mi(j)*rho;  
end
return
