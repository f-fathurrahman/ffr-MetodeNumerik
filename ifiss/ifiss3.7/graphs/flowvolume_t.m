function [uxref,total,upff] = flowvolume_t(xy,xns,xref,mkr,ttk,uxmax,uxmin)
%FLOWVOLUME_T plots flow solution on horizontal pipe section 
%   [uxref,total] = flowvolume_t(xy,xns,xref,fig,mkr,ttk,uxmax);
%   input
%          xy         velocity nodal coordinate vector 
%          xns        flow solution vector
%          xref       x-location of grid line  
%          mkr        plotting character
%          ttk        current time
%          uxmax      maximum x-velocity
%          uxmin      minimum x-velocity  
%   output
%          uxref      x-section horizontal flow vector 
%          uyref      x-section vertical flow vector 
%          upff       figure handle
%
%   IFISS function: DJS; 1 May 2012.
% Copyright (c) 2012 D.J. Silvester, H.C. Elman, A. Ramage 
% Modified  (c) 2022 M.D. Mihajlovic
nvtx=length(xy(:,1)); 
dd=(xy(:,1)-xref);
[m,ind]=min(abs(dd));
xref=xy(ind,1);
%fprintf('\nX-section analysis | x = %8.4e \n',xref)
kk=find(xy(:,1)==xref);
uxref=xns(kk);
yref=xy(kk,2);
npts=length(uxref);
upff=figure('visible','off');
plot(yref,uxref,mkr); 
axis([-1 1 1.05*uxmin 1.05*uxmax]);
title(['ux at x=' num2str(xref) ', t=' num2str(ttk)],'FontSize',10);
%%
%% compute volume of flow using appropriate quadrature
hy=diff(yref);
%% Simpson's rule
ww=ones(npts,1); ww(2:2:npts-1)=4*hy(1:2:npts-1); 
ww(1)=hy(1); ww(npts)=hy(npts-1);
ww(3:2:npts-2)=hy(2:2:npts-2)+hy(3:2:npts-2);
total=(1/3)*sum(ww.*uxref);
%fprintf('X-section flow volume is %8.4e \n',total)
return
