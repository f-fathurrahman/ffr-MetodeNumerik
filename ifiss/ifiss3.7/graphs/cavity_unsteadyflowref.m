function[fr,frx,fry,tfr]=cavity_unsteadyflowref...
  (qmethod,ev,usol,Tsol,tt,Av,BBy,BBx,Qv,xyv,x,y,bnd_d,npx,npy,net,xref,yref,hty,...
    i,krestart,vfreq,uvid,tvid,upvid,upxvid,upyvid,uxmax,uymax,fr,frx,fry,tfr)
%CAVITY_UNSTEADYFLOWREF plots cavity flow data at selected snapshot times 
%  [fr,frx,fry,tfr]=cavity_unsteadyflowref...
%   (qmethod,ev,usol,Tsol,tt,Av,BBy,BBx,Qv,xyv,x,y,bnd_d,npx,npy,net,xref,yref,hty,...
%    i,krestart,vfreq,uvid,tvid,upvid,upxvid,upyvid,uxmax,uymax,fr,frx,fry,tfr);
%   input
%          qmethod    mixed method 
%          mv2        Q2/Q1 element mapping matrix
%          usol       velocity solution (x+y)
%          Tsol       temperature solution
%          tt         discrete solution times
%          Av         vector diffusion matrix
%          BBy        velocity y-derivative matrix    
%          BBx        velocity x-derivative matrix    
%          Qv         velocity mass matrix
%          xyv        velocity nodal coordinate vector    
%          x          vector of x-axis interpolation points
%          y          vector of y-axis interpolation points
%          bnd_d      nodes on Dirichlet boundary
%          npx        number of elements in x direction
%          npy        number of elements in y direction
%          net        nodes ordering direction
%          xref       reference x grid line
%          yref       reference y grid line
%          hty        cavity problem type
%          i          segment number
%          krestart   number of segments
%          vfreq      plotting frequency
%          uvid       velocity video handler
%          tvid       temperature vodeo handler
%          upvid      velocity profile video handler
%          upxvid
%          upyvid
%          uxmax      maximum x velocity in the time segment
%          uymax      maximum y velocity in the time segment
%          fr         flow rates through the profile
%          frx        x-component flow rate
%          fry        y-component flow rate
%          tfr        discrete times of flow rates
%   output
%          fr         flow rates through the profile
%          tfr        discrete times of flow rates
%          frx        x-component flow rate
%          fry        y-component flow rate
% calls function xxstreambc.m to set boundary values
% Modified (c) 2022 M.D. Mihajlovic  
load box_bouss_plot.mat
warning off
%
nvtx=length(xyv); nu=2*nvtx;                  %  number of velocity dof
[LQv,UQv]=lu(Qv(1:nvtx,1:nvtx));              %  LU facrorisarion of vel. mass matrix
Asv=Av(1:nvtx,1:nvtx);                        %  scalar velocity Laplacian
fzero=zeros(nvtx,1);                          %  zero right-hand side
[Abc,fzero]=streambc(Asv,fzero,xyv,bnd_d);    %  impose zero BCs
[LA,UA]=lu(Abc);     %  LU factorisation of streamline Laplacian
%  loop over the snapshots
for k=1:vfreq:length(tt)
   ttk=tt(k);                                 %  selecting the discrete plotting times
%% Velocity streamlines
   u=usol(:,k);                               %  velocity at time tk
   ux=u(1:nvtx); uy=u(nvtx+1:nu);             %  velocity x and y component at tk
   utotal=sqrt(ux.*ux+uy.*uy);                %  2-norm of the velocity
   fsv=[BBy,-BBx]*u;                          %  right-hand side
   [fsv]=xxstreambc(Asv,fsv,xyv,bnd_d,ttk);   %  impose dynamic BC
   phi=UA\(LA\fsv);                           %  the stream function
   vff=figure('visible','off');
   copyobj(box_dom,vff);                      %  plot loop domain boundaries
   switch(net)
      case(1)    %  x direction lexicographical ordering 
         xx=reshape(xyv(:,1),2*npx+1,2*npy+1);      %  x-coordinates in grid format
         yy=reshape(xyv(:,2),2*npx+1,2*npy+1);      %  y-coordinates in grid format
      case(2)    %  y direction lexicographical ordering
         xx=reshape(xyv(:,1),2*npy+1,2*npx+1);      %  x-coordinates in grid format
         yy=reshape(xyv(:,2),2*npy+1,2*npx+1);      %  y-coordinates in grid format
   end
   xysol=griddata(xyv(:,1),xyv(:,2),phi,xx,yy); %  format the solution for plotting
   maxphi=max(max(xysol));                    %  max(phi)
   minphi=min(min(xysol));                    %  min(phi)
   contour(xx,yy,xysol,48);                   %  plot velocity streamlines
   axis([min(xyv(:,1))-0.2 max(xyv(:,1))+0.2 min(xyv(:,2))-0.2 max(xyv(:,2))+0.2]);
   title(sprintf('Velocity streamlines, t=%s',num2str(ttk)));
   vh=gca;
   framev=getframe(vff);
   writeVideo(uvid,framev);
   if(k>length(tt)-vfreq & k<=length(tt) & i==krestart)  %  save the plot at the end of the time interval
      vfig=figure(200);
      vgca=copyobj(vh,vfig); 
   else
      close(vff);
   end
%% Isotherms
   T=Tsol(:,k);                                %  teerature at time tk
   tff=figure('visible','off');
   copyobj(box_dom,tff);                       %  plot loop domain boundaries
   xysol=griddata(xyv(:,1),xyv(:,2),T,xx,yy);  %  format the solution for plotting
   contour(xx,yy,xysol,48);                    %  plot the isotherms
   axis([min(xyv(:,1))-0.2 max(xyv(:,1))+0.2 min(xyv(:,2))-0.2 max(xyv(:,2))+0.2]);
   title(sprintf('Temperature, t=%s',num2str(ttk)));
   th=gca;
   framet=getframe(tff);
   writeVideo(tvid,framet);
   if(k>length(tt)-vfreq & k<=length(tt) & i==krestart)  %  save the plot at the end of the interval
      tfig=figure(600);
      tgca=copyobj(th,tfig);
   else
      close(tff);
   end
 %% Velocity profile at specified lines
   switch(hty)   %  cavity problem type
      case(1)    %  Rayleigh-Benard
         [uxref,uyref,total,upff]=flowvolume_c(3,hty,0,xyv,u,xref,yref,'.-',ttk,uxmax,uymax);
         frameup=getframe(upff);
         writeVideo(upvid,frameup);
         close(upff);      
         fr=[fr,abs(total)];
         tfr=[tfr,ttk];
      case(2)    %  laterally heated cavity
         [uxref,uyref,total,uxpff]=flowvolume_c(3,hty,1,xyv,u,xref,yref,'.-',ttk,uxmax,uymax);
         frameup=getframe(uxpff);
         writeVideo(upxvid,frameup);
         close(uxpff);      
         frx=[frx,abs(total)];
         tfr=[tfr,ttk];
         [uxref,uyref,total,uypff]=flowvolume_c(3,hty,2,xyv,u,xref,yref,'.-',ttk,uxmax,uymax);
         frameup=getframe(uypff);
         writeVideo(upyvid,frameup);
         close(uypff);      
         fry=[fry,abs(total)];
         tfr=[tfr,ttk];
   end
end
return;
