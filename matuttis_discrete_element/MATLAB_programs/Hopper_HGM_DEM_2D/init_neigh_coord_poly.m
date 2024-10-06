function [a_r,b_r,c_r]=...
    init_neigh_coord_poly(partr, s)

n_part=length(s);
a_r=zeros(1,2*n_part+2);
b_r=zeros(1,2*n_part+2);
c_r=zeros(1,2*n_part+2);

% Sentinel with +-r_extrem to avoid unecessary ifs in later computation 
a_r(1)=-realmax; 
b_r(1)=0;
c_r(1)=0;
for i=1:n_part
  a_r(2*i)=min(partr(1:s(i),i));
  b_r(2*i)=-1;
  c_r(2*i)=i;
  a_r(2*i+1)=max(partr(1:s(i),i));
  b_r(2*i+1)=+1;
  c_r(2*i+1)=i;
end

a_r(2*n_part+2)=realmax;
b_r(2*n_part+2)=0;
c_r(2*n_part+2)=0;

return
