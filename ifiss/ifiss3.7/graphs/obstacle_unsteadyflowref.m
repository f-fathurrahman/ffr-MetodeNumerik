function obstacle_unsteadyflowref(qmethod,sol,tt,By,Bx,A,xy,x,y,...
						bound,bnde,xref,i,krestart,vfreq,uvid,upvid,uxmax,uxmin)
%OBSTACLE_UNSTEADYFLOWREF generates obstacle flow movie
%   obstacle_unsteadyflowref(qmethod,U,time,By,Bx,A,xy,x,y,bound,bnde,...
%                            xref,i,krestart,vfreq,uvid,upvid,uxmax,uxmin);
%   input
%          qmethod    mixed method 
%          U          velocity solution vector
%          time       solution time vector
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
% pressure solution is assumed to be essentially zero at outflow 
% so streamfunction satisfies zero Neumann condition there
% calls function xxstreambc.m to set boundary values
%   IFISS function: DJS; 3 August 2023.
% Copyright (c) 2009 D.J. Silvester, H.C. Elman, A. Ramage 
% Modified (c) 2023 M.D. Mihajlovic
load obst_ns_plot.mat
L=max(x); 
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
   u=sol(:,k);
   f=[By,-Bx]*u;
   [fsv]=xxstreambc(Asv,f,xy,bound,ttk);
   phi=UA\(LA\fsv);
   vff=figure('visible','off');
   copyobj(obst_dom,vff);
% interpolate to a cartesian product mesh
   [X,Y]=meshgrid(x,y);
% plot stream function
   colormap jet
   xysol = griddata(xy(:,1),xy(:,2),phi,X,Y);
   maxphi=max(max(xysol)); minphi=min(min(xysol));
%  if size(obs,1)~=0, [II,JJ] = findobsXY(obs,X,Y,bndxy); xysol(II,JJ)=nan; end
   vneg=[minphi:-minphi/12:0];
   vpos=[maxphi/12:maxphi/12:11*maxphi/12];
   vpospos=[maxphi/48: maxphi/48:maxphi/12];
   vnegneg=[minphi/12:-minphi/48:-minphi/48];
   contour(X,Y,xysol,[vneg,vpos,vpospos,vnegneg]); 
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
%%
return
