function [lo_bound,up_bound]=update_boundbox(s,r)
% PURPOSE: Update the bounding-boxes according to
%  the coordinates of the corners r of the current
%  timestep; the vector s contains the number of 
%  corners of each polygon
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014

n_part=length(s);
lo_bound=zeros(1,n_part);
up_bound=zeros(1,n_part);

for i = 1:n_part
   lo_bound(i) = 1000*min(r(1:s(i),i));
   up_bound(i) = 1000*max(r(1:s(i),i));
end

return
