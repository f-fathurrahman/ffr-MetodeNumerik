function [a]=update_a_list(lo_bound,up_bound,a,b,c)
% PURPOSE: Replace the bounding-box coordinates
%   from the previous timestep with the ones
%   in the current timestep in a_x,a_y
% USAGE: Called before the reordering in the 
%   neighborhood-algorithm
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014

lena=length(a);
for j=2:lena-1;
  i=c(j);
  if (b(j)<0)
    a(j)=lo_bound(i);
  else
    a(j)=up_bound(i);
  end
end

return
