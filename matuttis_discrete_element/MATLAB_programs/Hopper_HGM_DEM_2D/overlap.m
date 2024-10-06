function [overlap_geometry]=...
   overlap(p0x,p0y,pc0x,pc0y,p1x,p1y,pc1x,pc1y)
% PURPOSE: Compute the intersection geometry (overlap area, force
%  point, contact normal) of two overlapping polygonal
%  particles 
% ALGORITHM: Two nested loops over the edges
%  of the two polygons to computer the intersection
%  points; Additionally to intersecting edges,
%  the case of one edge touching the other has
%  also been treated
% CAVEAT:
%  - Algorithm will fail if polygons are  
%  - "Penetrating contacts" (sec. 7.3.6)
%     have not been dealt with
%  - The contact normal is not computed as in
%    eq. (7.9), but it is the normal to the 
%    line segment between the intersection points
%    of the outlines between the two polygons
% LITERATURE: sec. 7.3.2
% TODO: To speed up the interaction computation,
%  the information of the indices of the intersecting edges
%  from the previous timestep (p. 244) could be used to
%  to speed up the double loop computation
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014

% Normal 
nx=0;
ny=0;

maxNs =5;
maxNsp=maxNs*2;

sx=zeros(1,maxNs);
sy=zeros(1,maxNs);
np0s=zeros(1,maxNs);
np1s=zeros(1,maxNs);

np0=length(p0x);
np1=length(p1x);

% Correction in data structure of partcles.
% Build dummies.                            % Rev 1.
p0x=[p0x;p0x];
p0y=[p0y;p0y];
p1x=[p1x;p1x];
p1y=[p1y;p1y];

ns=0;
istart=1;
jstart=1;
iend  =np0;
jend  =np1;

% Loop over the edges of particle 0
for i=istart:iend
  r0x=p0x(i);
  r0y=p0y(i);
   
  u0x=p0x(i+1)-r0x;
  u0y=p0y(i+1)-r0y;
  p0minx=min( r0x, p0x(i+1) );
  p0miny=min( r0y, p0y(i+1) );
  p0maxx=max( r0x, p0x(i+1) );
  p0maxy=max( r0y, p0y(i+1) );
   
% Loop over the edges of particle 1   
  for j=jstart:jend
    r1x=p1x(j);
    if( r1x<p1x(j+1) )
      if( r1x>p0maxx || p1x(j+1)<p0minx )
        % goto_j_loopend=1; Rev 2.
        continue
      end
    else
      if(r1x<p0minx||p1x(j+1)>p0maxx)
        % goto_j_loopend=1; Rev 2.
        continue
      end
    end
    % end Rev 2.
      
    % if (goto_j_loopend~=1) Rev 2.
    r1y=p1y(j);
    if (r1y<p1y(j+1))
      if (r1y>p0maxy||p1y(j+1)<p0miny)
        % goto_j_loopend=1; Rev 2.
        continue
      end
    else
      if (r1y<p0miny||p1y(j+1)>p0maxy)
        % goto_j_loopend=1; Rev 2.
        continue
      end
    end
    % end Rev 2.
      
    % Possible overlap
    u1x=p1x(j+1)-r1x;
    u1y=p1y(j+1)-r1y;
    rx=r1x-r0x;
    ry=r1y-r0y;
    d0=u0y * u1x-u0x * u1y;
    if (d0>0)
      % First case
      d1=ry * u1x-rx * u1y;
      if (d1<=0||d1>d0)
        continue
      end
            
      d2=u0x*ry-u0y*rx;
      if (d2<=0||d2>d0)
        continue
      end       
      % Intersection point
      ns=ns+1;
      if (abs(d0)<1e-24)
        error('d0 soo small')
      end
      lambda=d1/d0;
    elseif(d0==0)
      if (rx*u0y~=ry*u0x)
        continue
      end     
% Sides on top of each other endpoint of the 0.
%  polygon selected as center of mass
      lambda=1.0d0;
      ns=ns+1;
    else
      % Second case
      d1=ry*u1x-rx*u1y;
      if (d1>=0||d1<d0)
        % goto_j_loopend=1; Rev 2.
        continue
      end
      d2=u0x*ry-u0y*rx;
      if (d2>=0||d2<d0)
        % goto_j_loopend=1; Rev 2.
         continue
      end     
      % Intersection point
      ns=ns+1;
      if (abs(d0)<1e-24)
        error ('d0 soo small')
      end
      lambda=d1/d0;
    end
    sx(ns)=r0x+lambda*u0x;
    sy(ns)=r0y+lambda*u0y;
    np0s(ns)=i;
    np1s(ns)=j;
    if (ns==1)
      inkoord(1)=i;
      inkoord(2)=j;
    end
    if (ns==2)
      inkoord(3)=i;
      inkoord(4)=j;
    end
  end  
end

sx=sx(1:ns);
sy=sy(1:ns);

spx=zeros(1,maxNsp);
spy=zeros(1,maxNsp);

nsp=0;

if ((ns<2)&&(ns>0)) % Rev. S.H. June 6, 2011
   fprintf('Only one intersection point???\n')
   fprintf('Impossible!!!\n')
   keyboard
elseif ( ns == 2 )
   % Two intersection points.
   % Good, now we create the overlap polygon.
   rx =sx(2)-sx(1);
   ry =sy(2)-sy(1);
   tmp=(pc1x-pc0x) * ry-(pc1y-pc0y) * rx;
   
   % First point.
   nsp=nsp+1;
   spx(nsp)=sx(1);
   spy(nsp)=sy(1);
   
 % TMP is computed (cross product) in order to get
 % correct orientation of the overlap area.
   if (tmp>0.d0)    
      % For polygon 0.
      ndif=np0s(2)-np0s(1)+1;
      if ( ndif < 0 )
         nadd=np0;
      else
         nadd=0;
      end
      
      for i=np0s(1)+1: np0s(2)+nadd
         nsp=nsp+1;
         spx(nsp)=p0x(i)-sx(1);
         spy(nsp)=p0y(i)-sy(1);
      end
      
      nsp=nsp+1;
      spx(nsp)=rx;
      spy(nsp)=ry;
      
      % For polygon 1.
      ndif=np1s(1)-np1s(2)+1;
      if (ndif<0)
         nadd=np1;
      else
         nadd=0;
      end
      
      for i=np1s(2)+1:np1s(1)+nadd
         nsp=nsp+1;
         spx(nsp)=p1x(i)-sx(1);
         spy(nsp)=p1y(i)-sy(1);
      end
      
      nx=sy(1)-sy(2);
      ny=sx(2)-sx(1);     
   else
      % For polygon 1.
      ndif=np1s(2)-np1s(1)+1;
      if (ndif<0)
         nadd=np1;
      else
         nadd=0;
      end
      
      for i=np1s(1)+1:np1s(2)+nadd
         nsp=nsp+1;
         spx(nsp)=p1x(i)-sx(1);
         spy(nsp)=p1y(i)-sy(1);
      end
      
      nsp=nsp+1;
      spx(nsp)=rx;
      spy(nsp)=ry;
      
      % For polygon 0.
      ndif=np0s(1)-np0s(2)+1;
      if (ndif<0)
         nadd=np0;
      else
         nadd=0;
      end
      
      for i=np0s(2)+1:np0s(1)+nadd
         nsp=nsp+1;
         spx(nsp)=p0x(i)-sx(1);
         spy(nsp)=p0y(i)-sy(1);
      end
      
      nx=sy(2)-sy(1);
      ny=sx(1)-sx(2);
   end
   
   if ((nx*nx+ny*ny)==0.d0)
      ns=0;
   end
   
elseif (ns>2)
  disp('There are more than 2 contact points!')
  error('Implement the algorithm as in sec. 7.2.6')
end

spx=spx(1:nsp);
spy=spy(1:nsp);

% Generate the geometry of the overlapped area.
if (nsp>0)
  overlap_area_x=zeros(1,nsp);
  overlap_area_y=zeros(1,nsp);
  overlap_area_x(1)=spx(1);
  overlap_area_y(1)=spy(1);
  for i=2:nsp
    overlap_area_x(i)=spx(1)+spx(i);
    overlap_area_y(i)=spy(1)+spy(i);
  end
  
  [scx, scy, area]=centerOfMass(overlap_area_x, overlap_area_y);
 
elseif (ns>2)
  [normalDir,forcePoint,sumArea]=weighted_normal(...
  sxCut0,syCut0,aCut0,sxCut1,syCut1,aCut1,...
  pc0x,pc0y,cxCut0,cyCut0,cxCut1,cyCut1);
  area=sumArea;
  scx =forcePoint(1);
  scy =forcePoint(2);
  nx  =normalDir(1);
  ny  =normalDir(2);
else   
  area=0;
  scx =0;
  scy =0;
  nx  =0;
  ny  =0;
end

overlap_geometry(1)=area;
overlap_geometry(2)=scx;
overlap_geometry(3)=scy;
overlap_geometry(4)=nx;
overlap_geometry(5)=ny;

return


