function step_unsteadytempref(sol,tt,xy,x,y,vfreq,i,krestart,tvid)
%STEP_UNSTEADYTEMPREF plots temperature data at snapshot times
%   step_unsteadytempref(sol,tt,xy,x,y,vfreq,i,tvid)
%   input
%          sol        temperature solution matrix (ordered by columns)
%          tt         solution time vector  
%          xy         velocity nodal coordinate vector   
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          vfreq      plotting fequency
%          i          the segment number
%          krestart   the number of segments
%          tvid       temperature video handle
%
%   IFISS function: DJS; 27 July 2015
% Copyright (c) 2011 D.J. Silvester, H.C. Elman, A. Ramage 
% Modified  (c) 2022 M.D. Mihajlovic
load step_bouss_plot.mat
L=max(x);        %  domain length
outflow=max(x);
warning off;
% interpolate to a cartesian product mesh
[X,Y]=meshgrid(x,y);
umax=max(sol(:,end)); umin=min(sol(:,1));
for k=1:vfreq:length(tt)
   ttk=tt(k);   %  discrete plotting time 
   xysol = griddata(xy(:,1),xy(:,2),sol(:,k),X,Y);
   [II,JJ]=find(X<0 & Y<0); xysol(II,JJ)=nan;
   tff=figure('visible','off');
   copyobj(step_dom,tff);                     %  plot loop domain boundaries
   colormap jet
   contour(X,Y,xysol,30), axis('equal')
   axis([-1,outflow,-1,1]), stepx, axis('off')
   title(['Temperature: t=',num2str(ttk,'%5.2f')],'FontSize',12);
   th=gca;
   framet=getframe(tff);
   writeVideo(tvid,framet);
   if(k>length(tt)-vfreq & k<=length(tt) & i==krestart)     %  save the plot at the end of the time interval
      tfig=figure(800);
      tgca=copyobj(th,tfig);
   else
      close(tff);
   end
end
return
