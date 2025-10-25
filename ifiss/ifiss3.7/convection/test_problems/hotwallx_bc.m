
function bc = specific_bc(xbd,ybd)
%hotwallx_bc   Reference problem 3.4 regularized boundary condition 
%   bc = specific_bc(xbd,ybd);
%   input
%          xbd          x boundary coordinate vector
%          ybd          y boundary coordinate vector 
%
%   specifies regularised hot wall /Morton pp.10/
%   IFISS function: DJS; 24 May 2019.
% Copyright (c) 2005 D.J. Silvester, H.C. Elman, A. Ramage 
nobd=length(xbd);
bc=zeros(size(xbd));
k=find(xbd==1);
bc(k)=(1-ybd(k).*ybd(k)).*(1+ybd(k).*ybd(k));
return
