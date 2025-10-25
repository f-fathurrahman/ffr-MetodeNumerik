function channel_unsteadyflowref(qmethod,U,tt,By,Bx,A,xy,x,y,...
						bound,bnde,xref,i,krestart,vfreq,uvid,upvid,uxmax,uxmin)
%CHANNEL_UNSTEADYFLOWREF generates channel flow movie
%   channel_unsteadyflowref(qmethod,U,tt,By,Bx,A,xy,x,y,bound,bnde,...
%                           xref,i,krestart,vfreq,uvid,upvid,uxmax,uxmin);
%   input
%          qmethod    mixed method 
%          U          velocity solution vector
%          tt         discrete solution times vector
%          By         velocity  y-derivative matrix    
%          Bx         velocity x-derivative matrix    
%          A          vector diffusion matrix
%          xy         velocity nodal coordinate vector   
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          bound      boundary vertex vector
%          bnde       boundary edges vector
%          xref       reference x grid line
%          i          segment number
%          krestart   number of segments
%          vfreq      plotting frequency
%          uvid       velocity video handler
%          upvid      velocity profile video handler
%          uxmax      maximum x velocity in the time segment  
%          uxmin      minimum x velocity in the time segment  
%          
%
% pressure solution is assumed to be essentially zero at outflow 
% so streamfunction satisfies zero Neumann condition there
% calls function xxstreambc.m to set boundary values
%   IFISS function: DJS; 3 August 2023.
% Copyright (c) 2009 D.J. Silvester, H.C. Elman, A. Ramage 
% Modified (c) 2023 M.D. Mihajlovic
load channel_plot.mat
nvtx=length(xy); nu=2*nvtx; 
Asv=A(1:nvtx,1:nvtx); fzero=zeros(nvtx,1);
%
%% apply boundary conditions to the diffusion matrix
[Abc,fzero]=streambc(Asv,fzero,xy,bound);
[LA,UA]= lu(Abc); 
%
% ------------------ loop over snapshots
for k=1:vfreq:length(tt)
   ttk=tt(k); 
%% Velocity streamlines  
   u=U(:,k);
   f=[By,-Bx]*u;
   [fsv]=xxstreambc(Asv,f,xy,bound,ttk);
   phi=UA\(LA\fsv);
%  fprintf('  %4i   %7.3f  %10.5f   %9.3e\n', kk, ttk, min(phi),max(phi));
   vff=figure('visible','off');
   copyobj(channel_dom,vff);
% interpolate to a cartesian product mesh
   [X,Y]=meshgrid(x,y);
% plot stream function
   colormap jet
   xysol = griddata(xy(:,1),xy(:,2),phi,X,Y);        %  format the solution for plotting
   maxphi=max(max(xysol)); minphi=min(min(xysol));   %  max(phi) and min(phi)
   vneg=[minphi:-minphi/48:0];
   vpos=[maxphi/48:maxphi/48:47*maxphi/48];
   contour(X,Y,xysol,[vneg,vpos]);
   %contour(X,Y,xysol,48);
%  axis([min(x)-0.2 max(x)+0.2 min(y)-0.2 max(y)+0.2]); 
   title(['Velocity streamlines: time = ',num2str(ttk,'%6.3f')],'FontSize',10),
   vh=gca;
   framev=getframe(vff);
   writeVideo(uvid,framev);
   if(k>length(tt)-vfreq & k<=length(tt) & i==krestart)  %  save the plot at the end of the time interval
      vfig=figure(200);
      vgca=copyobj(vh,vfig); 
   else
      close(vff);
   end
%% Velocity profile at specified lines
   [uxref,total,upff]=flowvolume_t(xy,u,xref,'.-',ttk,uxmax,uxmin);
   frameup=getframe(upff);
   writeVideo(upvid,frameup);
   close(upff);
end
% ------------------ end loop over snapshots

return
