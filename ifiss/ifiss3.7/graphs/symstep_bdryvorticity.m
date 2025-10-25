function [label,it,ib] = symstep_bdryvorticity(domain,qmethod,xy,bound,w,fig,mkr)
%SYMSTEP_BDRYVORTICITY plots vorticity solution on symmetric step boundary
%   [label,it,ib] = symstep_bdryvorticity(domain,qmethod,xy,bound,w,19,'or-');
%   input
%          domain     problem domain 
%          qmethod    mixed method 
%          xy         velocity nodal coordinate vector 
%          bound      boundary node index 
%          w          vorticity solution vector
%          fig        figure reference 
%          mkr        plotting character
%   output
%          label      1/0 symmetric/nonsymmetric label
%             it      top boundary reattachment index
%             ib      bottom boundary reattachment index
%
%   IFISS function: DJS; 17 April 2023
% Copyright (c) 2012 D.J. Silvester, H.C. Elman, A. Ramage
nvtx=length(xy(:,1)); nbound=length(bound);
fprintf('\n vorticity boundary analysis\n')
if domain~=3, error('function is only defined for symmetric step'), end

% ------- lower boundary
kk=find(xy(bound,2)==-1);
nkk=length(kk);wbound=w(bound(kk));xbound=xy(bound(kk),1);
figure(fig)
subplot(121)
plot(xbound,wbound,'-bo'), axis('square'),
xlabel('x');
title('lower boundary vorticity distribution','FontSize',11)
ib=find(wbound(5:end)>0,1)
[wmin,kmin]=min(wbound);
fprintf('minimum vorticity is %8.4e \n',wmin)
fprintf('           when x is  %8.4e \n',xbound(kmin))
%
% ------- upper boundary
kk=find(xy(bound,2)==1);
nkb=length(kk);wboundb=w(bound(kk));xbound=xy(bound(kk),1);
subplot(122)
plot(xbound,wboundb,mkr), axis('square'),
xlabel('x');
title('upper boundary vorticity distribution','FontSize',11)
it=find(wboundb(5:end)<0,1)
[wmax,kmax]=max(wboundb);
fprintf('maximum vorticity is %8.4e \n',wmax)
fprintf('           when x is  %8.4e \n',xbound(kmax))
if abs(ib-it)<2 % indices within 1 of each other
   label = 0;
   fprintf('flow solution is SYMMETRIC \n')
else
   label = 1;
   fprintf('flow solution is NONSYMMETRIC \n')
end
return
