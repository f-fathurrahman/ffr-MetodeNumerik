function [uxref,uyref,total,upff] = flowvolume_c(qmethod,hty,cmp,xy,xns,xref,yref,mkr,ttk,uxmax,uymax)
%FLOWVOLUME_C plots flow solution on horizontal or vertical section of a cavity
%   [uxref,uyref,total,upff] = flowvolume_c(qmethod,hty,cmp,xy,xns,xref,yref,mkr,ttk,uxmax,uymax);
%   input
%          qmethod    mixed method
%          hty        cavity type problem
%          cmp        velocity component
%          xy         velocity nodal coordinate vector 
%          xns        flow solution vector
%          xref       x location of a grid line  
%          yref       y location of a grid line
%          mkr        plotting character
%          ttk        current time
%          uxmax      maximum x-velocity
%          uymax      maximum y-velocity
%   output
%          uxref      x-component velocity vector 
%          uyref      y-component velocity vector 
%          total      flow rate at the intersection 
%          upff       figure handle
%
%   IFISS function: DJS; 1 May 2012.
% Copyright (c) 2012 D.J. Silvester, H.C. Elman, A. Ramage 
% Modified  (c) 2022 M.D. Mihajlovic
upff=figure('visible','off');
nvtx=length(xy(:,1));
switch(hty)  %  cavity problem switch
   case(1)   %  Rayleigh-Benard
      kk=find(xy(:,2)==yref);
      uxref=xns(kk);
      uyref=xns(nvtx+kk); 
      npts=length(uyref);
      xya=xy(kk,1);
      plot(xya,uyref,mkr); 
      axis([min(xya) max(xya) -uymax uymax]);
      title(['uy at y=' num2str(yref) ', t=' num2str(ttk)],'FontSize',10);
   case(2)   %  laterally heated cavity
      switch(cmp)   %  velocity component switch
          case(1)   %  x-component 
             kk=find(xy(:,1)==xref);
             uxref=xns(kk);
             uyref=xns(nvtx+kk); 
             npts=length(uxref);
             xya=xy(kk,2);
             plot(xya,uxref,mkr); 
             axis([min(xya) max(xya) -uxmax uxmax]);
             title(['ux at x=' num2str(xref) ', t=' num2str(ttk)],'FontSize',10); 
         case(2)    %  y-component
             kk=find(xy(:,2)==yref);
             uxref=xns(kk);
             uyref=xns(nvtx+kk); 
             npts=length(uyref);
             xya=xy(kk,1);
             plot(xya,uyref,mkr); 
             axis([min(xya) max(xya) -uymax uymax]);
             title(['uy at y=' num2str(yref) ', t=' num2str(ttk)],'FontSize',10);
      end
end
%%
%% compute volume of flow using appropriate quadrature
gref=xya;
hy=diff(gref);
%% Simpson's rule
ww=ones(npts,1); ww(2:2:npts-1)=4*hy(1:2:npts-1); 
ww(1)=hy(1); ww(npts)=hy(npts-1);
ww(3:2:npts-2)=hy(2:2:npts-2)+hy(3:2:npts-2);
if(hty==1)
   total=(1/3)*sum(ww.*uyref);
else
   switch(cmp)
      case(1)  %  x-velocity   
         total=(1/3)*sum(ww.*uxref);
      case(2) %  y-velocity
         total=(1/3)*sum(ww.*uyref);
   end
end
return
