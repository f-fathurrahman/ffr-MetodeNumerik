function square_unsteadyflowref(qmethod,U,tt,ev,By,Bx,A,G,xy,xyp,x,y,bound,...
                                 xref,yref,i,krestart,vfreq,uvid,vvid,uxpvid,uypvid,uxmax,uxmin,uymax,uymin)
%SQUARE_UNSTEADYFLOWREF  generates cavity flow movie
%   square_unsteadyflowref(qmethod,U,tt,ev,By,Bx,A,G,xy,xyp,x,y,bound,...
%                          xref,yref,i,krestart,vfreq,uvid,vvid,uxpvid,uypvid,uxmax,uxmin,uymax,uymin);
%   input
%          qmethod    mixed method 
%          U          flow solution vector
%          tt         snapshot time vector
%          ev         mapping
%          By         velocity  y-derivative matrix    
%          Bx         velocity x-derivative matrix   
%          A          vector diffusion matrix
%          G          veclocity mass matrix
%          xy         velocity nodal coordinate vector  
%          xyp        pressure nodal coordinate vector  
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          bound      boundary vertex vector
%          xref       reference x grid line
%          yref       reference y grid line
%          i          segment number
%          krestart   number of segments
%          vfreq      plotting frequency
%          uvid       velocity video handler
%          vvid       vorticity video handler
%          uxpvid     x velocity profile video handler
%          uypvid     y velocity profile video handler
%          uxmax      maximum x velocity in the time segment  
%          uxmin      minimum x velocity in the time segment  
%          uymax      maximum y velocity in the time segment
%          uymin      minimum y velocity in the time segment
%
% calls function xxstreambc.m to set boundary values
%   IFISS function: DJS; 3 August 2023.
% Copyright (c) 2012 D.J. Silvester, H.C. Elman, A. Ramage 
% Modified (c) 2023 M.D. Mihajlovic
load cavity_ns_plot.mat
nvtx=length(xy); nu=2*nvtx; np=length(xyp);
[LG,UG]= lu(G(1:nvtx,1:nvtx)); 
Asv=A(1:nvtx,1:nvtx); fzero=zeros(nvtx,1);
%
%% apply boundary conditions to diffusion matrix
[Abc,fzero]=streambc(Asv,fzero,xy,bound);
[LA,UA]=lu(Abc); 
%
% ------------------ loop over snapshots
for k=1:vfreq:length(tt)
   ttk=tt(k);
%% Velocity streamlines
   u=U(:,k);
   ux=u(1:nvtx); uy=u(nvtx+1:nu);  utotal=sqrt(ux.*ux+ uy.*uy);
   fsv=-[By,-Bx]*u;
   omega=UG\(LG\fsv);
   f=[By,-Bx]*u;
   [fsv]=xxstreambc(Asv,f,xy,bound,ttk);   
   phi=UA\(LA\fsv);  
%% Vorticity   
   if qmethod > 1, wev = vorticity_q2(xy,ev,omega,0);
   else, wev = vorticity_q1(xy,ev,omega,0); end
%
   uff=figure('visible','off');
   copyobj(cavity_dom,uff);
%  interpolate to a cartesian product mesh
   [X,Y]=meshgrid(x,y);
%  plot stream function
   colormap jet
   xysol = griddata(xy(:,1),xy(:,2),phi,X,Y);        %  format the solution for plotting  
   maxphi=max(max(xysol)); minphi=min(min(xysol));   %  max(phi) and min(phi)
   vneg=[minphi:-minphi/24:0];
   vpos=[maxphi/6:maxphi/6:maxphi];
   vpospos=[0: maxphi/48:maxphi/12];
   contour(X,Y,xysol,[vneg,vpos,vpospos])
   title(['Streamlines: time = ',num2str(ttk,'%5.2f')],'FontSize',12),
   uh=gca;
   framev=getframe(uff);
   writeVideo(uvid,framev);
   if(k>length(tt)-vfreq & k<=length(tt) & i==krestart)  %  save the plot at the end of the time interval
      ufig=figure(200);
      ugca=copyobj(uh,ufig); 
   else
      close(uff);
   end
%  plot vorticity
   vff=figure('visible','off');
   copyobj(cavity_dom,vff);
   xysol = griddata(xy(:,1),xy(:,2),omega,X,Y);
   solheight = max(max(xysol))-min(min(xysol));
   contour(X,Y,xysol,48)
   title(['Vorticity: time = ',num2str(ttk,'%5.2f')],'FontSize',12),  
   vh=gca;
   framev=getframe(vff);
   writeVideo(vvid,framev);
   if(k>length(tt)-vfreq & k<=length(tt) & i==krestart)  %  save the plot at the end of the time interval
      vfig=figure(300);
      vgca=copyobj(vh,vfig); 
   else
      close(vff);
   end
%% Horizontal velocity profile
   [uxref,uyref,total,upff]=flowvolume_c(qmethod,2,1,xy,u,xref,yref,'.-',ttk,uxmax,uymax);
   frameup=getframe(upff);
   writeVideo(uxpvid,frameup);
   close(upff);
%% Vertical velocity profile
   [uxref,uyref,total,upff]=flowvolume_c(qmethod,2,2,xy,u,xref,yref,'.-',ttk,uxmax,uymax);
   frameup=getframe(upff);
   writeVideo(uypvid,frameup);
   close(upff);
end
% ------------------ end loop over snapshots
return
