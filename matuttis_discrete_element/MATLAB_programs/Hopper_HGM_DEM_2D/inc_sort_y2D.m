function [wrong_position_index_2,wrong_position_list_2,a,b,c]=... 
inc_sort_y2D(wrong_position_list_1,wrong_position_index_1,a,b,c)
% PURPOSE: Incrementally sort the bounding boxes by binary 
%   exchange (bubble-sort) output a list of possible candidate 
%   positions where reordering is still necessary
% USAGE: Can be used only for the y-direction, 
%   for the x-direction, additional conditions are to
%   be taken into account
% CAVEAT: 
%  - lo_bound_x_old AND up_bound_x_old must be
%   Initialized as 0 in first step !!!!!
%  - If at the inital timestep, the particles
%   have been initialiyed in unfavorable positions 
%   (particles with lower index above particles with
%   higher index), for n particles the algorithm may take
%   time proportional to n*n, see p. ????
% REVISION HISTORY: HG Matuttis, Wed Aug 18 13:48:16 MSZ 1999
%
global new_col_list
global new_col_list_len
global lo_bound_x;
global up_bound_x;
global lo_bound_y;
global up_bound_y;

global lo_bound_x_old;
global up_bound_x_old;

len_a=length(a);

%1. Binary exchange of wrongly ordered elements, saving of
%   indices where ordering may be necessary
following_position_index_1=0;
for wrong_list_index=1:wrong_position_index_1 % verified ! start at 1, not 2
   first_index=wrong_position_list_1(wrong_list_index);
   secon_index=first_index+1;
% exchange particles
   dummy=a(first_index);
   a(first_index)=a(secon_index);
   a(secon_index)=dummy;
   dummy=b(first_index);
   b(first_index)=b(secon_index);
   b(secon_index)=dummy;
   dummy=c(first_index);
   c(first_index)=c(secon_index);
   c(secon_index)=dummy;
% Set up collision list:
   if ( (c(first_index)~=c(secon_index))&&...                    % 1.
        (b(first_index)*b(secon_index)<0)&&...                   % 2.  
       (((lo_bound_x(c(first_index))<up_bound_x(c(secon_index)))&&...% 3.
         (lo_bound_x(c(secon_index))<up_bound_x(c(first_index)))))&&...         
       (((lo_bound_y(c(first_index))<up_bound_y(c(secon_index)))&&...% 4.
         (lo_bound_y(c(secon_index))<up_bound_y(c(first_index))))))
% 1. Avoid selfinteraction
% 2. Only if "correct" neighbors cross
% 3. Only if total siutation is like below:
% |--- i---|       lo_bound_y(i)< up_bound_y(j)&lo_bound_y(j)<up_bound_y(i)
%   |--- j---|             is the same as
%      |--- i---|  lo_bond(j) < up_bound_y(i) & lo_bound_y(i)<up_bound_y(j)
% Only inclusion in the collision list, if
% - the boundary boxes are the right ones AND
% - the coordinates have the right relative position
% -> WORKS ALSO FOR JUMPING PARTICLES, WITHOUT CONTINOUS MOTION
% Collision in 1D
  
     if ((lo_bound_x_old(c(first_index))<up_bound_x_old(c(secon_index)))&&...
         (lo_bound_x_old(c(secon_index))<up_bound_x_old(c(first_index))))     
       new_col_list_len=new_col_list_len+1;
       new_col_list(1,new_col_list_len)=min(c(first_index),c(secon_index));
       new_col_list(2,new_col_list_len)=max(c(first_index),c(secon_index));
     end
   end
% list with first_index-1,secon_index, give the indices of the left 
% particles in a comparison
% where errors can occur after the interchange
% No if-conditions like  (a(secon_index)>a(secon_index+1))
% because that would inhibit vectorization
% The lower index of particles which must be compared will be saved
   following_position_index_1=following_position_index_1+1; 
   following_position_list_1(following_position_index_1)=first_index-1;
   following_position_index_1=following_position_index_1+1; 
   following_position_list_1(following_position_index_1)=secon_index;
end

%2. Eliminate entries in the wrong order and
% redundant entries from the index list
reduced_following_position_index_1=0;
reduced_following_position_list_1=zeros(1,len_a);
% Avoid unnecessary special treatment of last entry
following_position_list_1(following_position_index_1+1)=1e200;
for i=1:following_position_index_1
  if ~(following_position_list_1(i)>=following_position_list_1(i+1))
    reduced_following_position_index_1=...
    reduced_following_position_index_1+1;
    reduced_following_position_list_1(reduced_following_position_index_1)=...
    following_position_list_1(i);
  end
end

%3. Set up list of particles which must be still reordered
wrong_position_index_2=0;
wrong_position_list_2=zeros(1,len_a);
for i=1:reduced_following_position_index_1
  j=reduced_following_position_list_1(i);
  if (a(j)>a(j+1))
    wrong_position_index_2=wrong_position_index_2+1;
    wrong_position_list_2(wrong_position_index_2)=j;
  end
end

return
