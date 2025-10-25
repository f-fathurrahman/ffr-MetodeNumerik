function step_unsteadyflowref(qmethod,mv2,U,tt,Av,BBy,BBx,Qv,xyv,xyp,x,y,bound,i,krestart,vfreq,xref,...
                                 uvid,vvid,upvid,uxmax,symm)
%STEP_UNSTEADYFLOWREF generates flow over step movie
% step_unsteadyflowref(2,mv2,U,soltime,Av,BBy,BBx,Qv,xyv,xyp,x,y,bound,i,vfreq,xref,uvid,vvid,uxmax,symm);
%   input
%          qmethod    mixed method 
%          ev         mv/ev  Q2/Q1 element mapping matrix
%          U          flow solution vector
%          tt         snapshot time vector
%          A          vector diffusion matrix
%          By         velocity  y-derivative matrix    
%          Bx         velocity x-derivative matrix    
%          G          veclocity mass matrix
%          xy         velocity nodal coordinate vector  
%          xyp        pressure nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          bound      boundary vertex vector
%          i          the segment number
%          krestart   the number of segments
%          vfreq      plotting frequency
%          xref       reference grid line
%          uvid       velocity streamlines video handle
%          vvid       vorticity video handle
%          upvid      velocity profile video handle
%          uxmax      maximum x-velocity in the segment
%          symm       symmetric domain 0/1 switch (optional)
%
% calls function xxstreambc.m to set boundary values
% calls function flowvolume_c.m to calculate the flow rate
%   IFISS function: DJS; 3 August 2023.
% Copyright (c) 2014 D.J. Silvester, H.C. Elman, A. Ramage
% Modified  (c) 2022 M.D. Mihajlovic
if nargin < 21, symm = 0; end
load step_bouss_plot.mat
L=max(x);    %  domain length
nvtx=length(xyv); nu=2*nvtx; np=length(xyp);
[LG,UG]=lu(Qv(1:nvtx,1:nvtx)); 
Asv=Av(1:nvtx,1:nvtx); fzero=zeros(nvtx,1);
[Abc,fzero]=streambc(Asv,fzero,xyv,bound);
[LA,UA]=lu(Abc); 
%
% ------------------ loop over snapshots
for k=1:vfreq:length(tt)
   ttk=tt(k);           %  selecting discrete plotting times
% compute derived quantites
   u=U(:,k);            %  velocity at time tk
   ux=u(1:nvtx);        %  x-velocity at tk
   uy=u(nvtx+1:nu);     %  y-velocity at tk
   utotal=sqrt(ux.*ux+ uy.*uy);   %  norm of the velocity
   fsv=-[BBy,-BBx]*u;   %  vorticity system right-hand side
   omega=UG\(LG\fsv);   %  vorticity
   f=[BBy,-BBx]*u;      %  streamline system righ-hand side
   [fsv]=xxstreambc(Asv,f,xyv,bound,ttk);  %  BCs 
   phi=UA\(LA\fsv);     %  streamline function
%  if(qmethod > 1)
%     wev = vorticity_q2(xyv,ev,omega,0);   %  Q2 approximation
%  else   
%     wev = vorticity_q1(xyp,ev,omega,0);   %  Q1 approximation
%  end
%
% interpolate to a cartesian product mesh
   [X,Y]=meshgrid(x,y);
%
% plot stream function
   uff=figure('visible','off');
   copyobj(step_dom,uff);                            %  plot loop domain boundaries
   ax = [min(x)-.2 max(x)+.2 min(y)-.2 max(y)+.2];   %  axis range
   xysol = griddata(xyv(:,1),xyv(:,2),phi,X,Y);      %  format the solution for plotting
   if(symm)    % symmetric step
      [II,JJ]=find(X<0 & Y<-0.5); xysol(II,JJ)=nan;
      [II,JJ]=find(X<0 & Y>0.5); xysol(II,JJ)=nan;
   else
      [II,JJ]=find(X<0 & Y<0); xysol(II,JJ)=nan;
   end
   maxphi=max(max(xysol)); minphi=min(min(xysol));   %  range for the streamfunction
   vneg=[minphi:-minphi/6:0];
   vpos=[maxphi/20:maxphi/20:19*maxphi/20];
   vpospos=[79*maxphi/80: maxphi/320:maxphi];
   if L<=5 %default domain
      contour(X,Y,xysol,[vneg,vpos]) 
      axis equal
   else
      contour(X,Y,xysol,[vneg,vpos,vpospos])
      %contour(X,Y,xysol,[vneg,vpospos])
      axis equal, axx=ax(1:4); axx(2)=min(30,axx(2));
      axis(axx(1:4)); 
   end
   title(['Velocity streamlines: t=',num2str(ttk,'%6.2f')],'FontSize',12), 
   if symm, stepsym, else, stepx, end, axis('off')
   uh=gca;
   frameu=getframe(uff);
   writeVideo(uvid,frameu);
   if(k>length(tt)-vfreq & k<=length(tt) & i==krestart)  %  save the plot at the end of the time interval
      ufig=figure(200);
      ugca=copyobj(uh,ufig);
   else
      close(uff);
   end
%
% plot vorticity
   vff=figure('visible','off');
   copyobj(step_dom,vff);                          %  plot loop domain boundaries
   xysol = griddata(xyv(:,1),xyv(:,2),omega,X,Y);  %  format the vorticity for plotting
   solheight = max(max(xysol))-min(min(xysol));    %  vorticity range
   ax5 = min(min(xysol))-.1*solheight;
   ax6 = max(max(xysol))+.1*solheight;
   ax = [min(x)-.2 max(x)+.2 min(y)-.2 max(y)+.2 ax5 ax6];   %  axis range
   if(symm)    % symmetric step
      [II,JJ]=find(X<0 & Y<-0.5); xysol(II,JJ)=nan;
      [II,JJ]=find(X<0 & Y>0.5); xysol(II,JJ)=nan;
   else
      [II,JJ]=find(X<0 & Y<0); xysol(II,JJ)=nan;
   end
   maxphi=max(max(xysol)); minphi=min(min(xysol));   %  range for vorticity
   vneg=[minphi:-minphi/12:0];
   vpos=[maxphi/20:maxphi/20:19*maxphi/20];
   vpospos=[79*maxphi/80: maxphi/320:maxphi];
   if L<=5 %default domain
      contour(X,Y,xysol,30)    
      %contour(X,Y,xysol,[vneg,vpos])
      axis equal
   else
      contour(X,Y,xysol,30)       
      %contour(X,Y,xysol,[vneg,vpos,vpospos])
      %contour(X,Y,xysol,[vneg,vpospos])
      axis equal, axx=ax(1:4); axx(2)=min(25,axx(2));
      axis(axx(1:4)); 
   end
   if symm, stepsym, else, stepx, end, axis('off'), axis('off')
   vh=gca;
   title(['Vorticity contours: t=',num2str(ttk,'%6.2f')],'FontSize',12);
   framev=getframe(vff);
   writeVideo(vvid,framev);
   if(k>length(tt)-vfreq & k<=length(tt) & i==krestart)   %  save the plot at the end of the time interval
      vfig=figure(400);
      ugca=copyobj(vh,vfig);
   else
      close(vff);
   end
%
%  plot velocity profile
   [uxref,uyref,total,uph]=flowvolume_c(3,2,1,xyv,u,xref,0,'.-',ttk,uxmax,0);
   frameup=getframe(uph);
   writeVideo(upvid,frameup);
   close(uph);
end
% ------------------ end loop over snapshots
%
return
