function [col_list,col_list_len]=...
  neighborhood_algorithm(polys,polyx,polyy)
%  PURPOSE: Neighborhood algorithm via "sort and sweep"
%    1. Re-compute the bounding boxes
%    2. Update the corresponding lists A_X, AY
%    3. Permute the lists A_X, ... A_Y ..., C_Y and 
%       register possible interacting particle pairs 
%       according to the overlap of bounding boxes
% INPUT:
%  polys - Number of corners in each polygon.
%  polyx - Corners of each polygons.
%  polyy
% OUTPUT:
%  col_list - List of possible contacting polygons. 
%  col_list_len - Number of possible contacting polygons.
% USAGE: The persistent variables in this function must be
%  cleared at the very beginning of the main program with
%   clear neighborhood_algorithm
% REVISION HISTORY:
%  11-May-2014 HG Matuttis Simplified interfaces
%  2013-Jun-21 S. H. Ng       Added into DEM-code.
%  1-2000      HG Matuttis Original
%
% ALGORITHM: Sort and Sweep
%  Vectorizable (no bank conflicts would arise on
%  a NEX-SX supercomputer)incremental sort bounding boxes:
%  x-sorting: take into account only new bounding boxes, like in 1D.
%  y-sorting: take into account old and new bounding boxes.
%
% CAVEAT: If the bounding boxes of some particles
%  overlap or even TOUCH during initialization, the 
%  neighborhood algorithm will not be able to detect
%  these interacting pairs
% TODO: Write a control routine for the particle
%  initialization so that overlaps during the 
%  initialization are detected to avoid the problem
%  in "CAVEAT"
% LITERATURE: sec. 7.5
% REVISION HISTORY: Shi Han Ng, HG Matuttis 15-May-2014

% Global variables for neighborhood algorithm.
persistent ISINITIALIZED
persistent A_X B_X C_X
persistent A_Y B_Y C_Y
persistent OLD_COL_LIST OLD_COL_LIST_LEN
global new_col_list
global new_col_list_len
global lo_bound_x;
global up_bound_x;
global lo_bound_y;
global up_bound_y;
global lo_bound_x_old;
global up_bound_x_old;

polyn=length(polys);

% Initialization of the neighborhood algorithm
% if the function is called for the first time
if (isempty(ISINITIALIZED) )
  OLD_COL_LIST_LEN=0;
  new_col_list=zeros(2,8*polyn);
  OLD_COL_LIST=zeros(2,8*polyn);
 
% Initialize the x- and y- coordinate.
  [A_X,B_X,C_X]=init_neigh_coord_poly(polyx,polys);
  lo_bound_x=zeros(1,polyn);
  up_bound_x=zeros(1,polyn);
  [A_Y,B_Y,C_Y]=init_neigh_coord_poly(polyy,polys);
  disp('Initialization of the neighborhood algorithm finished')
  ISINITIALIZED=1;
end

% Save old x-bounding boxes, only they are needed
lo_bound_x_old=lo_bound_x;
up_bound_x_old=up_bound_x;
% Update the bounding boxes.
[lo_bound_x,up_bound_x]=update_boundbox(polys,polyx);
[lo_bound_y,up_bound_y]=update_boundbox(polys,polyy);

% Update the lists A_X, A_X
[A_X]=update_a_list(lo_bound_x,up_bound_x,A_X,B_X,C_X);
[A_Y]=update_a_list(lo_bound_y,up_bound_y,A_Y,B_Y,C_Y);
   
new_col_list_len=0;
% Make a list of particles newly in contact,
% new_col_list, which is defined as
% global variable to avoid new initialization
[A_X,B_X,C_X]=start_search_x2D(A_X,B_X,C_X);
[A_Y,B_Y,C_Y]=start_search_y2D(A_Y,B_Y,C_Y);

% Append those elements to the new collision lists
% from the old collision lists where the bounding
% boxes still have an overlap
for i_list=1:OLD_COL_LIST_LEN
   i=OLD_COL_LIST(1,i_list);
   j=OLD_COL_LIST(2,i_list);
   if ((lo_bound_x(i)<up_bound_x(j))&&(lo_bound_x(j)<up_bound_x(i)))&&...
      ((lo_bound_y(i)<up_bound_y(j))&&(lo_bound_y(j)<up_bound_y(i)))
      new_col_list_len=new_col_list_len+1;
      new_col_list(1,new_col_list_len)=i;
      new_col_list(2,new_col_list_len)=j;
   end
end
OLD_COL_LIST=new_col_list;
OLD_COL_LIST_LEN=new_col_list_len;
col_list=OLD_COL_LIST;
col_list_len=OLD_COL_LIST_LEN;

return