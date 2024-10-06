function [a,b,c,search_depth,num_perm]=start_search_y2D(a,b,c)
%function [a,b,c,search_depth,num_perm]=start_search_y(a,b,c)
% HG Tue Dec 21 10:32:57 JST 1999
% PURPOSE: sort single array  using a vectorizable sorting algorithm
%   calls inc_sort_y
% USAGE: Prepare data so that inc_sort can be used
% ALGORITHM: Sorting by binary permutations of a list
% CAVEAT: The data in the input array aold must be confined
%   between extremal positions which act as "sentinel" 
%   so that unnecessary if-conditions for the end-of array-determination
%   can be omitted
% VARIABLES:
%   a: Array to be sorted
%    num_perm=number of permutations
% ALGORITHMS: Use two calls of inc_sort so that no
%   dummy copying of arrays is necessary
% ALGORITHM: 
% - Use two calls of inc_sort so that no
%  dummy copying of arrays is necessary
% OUT-COMMENTED LINES:
% - Allow to analyse the computational effort, which 
%  will vary with the physical system.
%  search_depth: Counts how often inc_sort_x2D is 
%  called.
%  num_perm: adds up the number of permuations
% - Control sums for the coordinates
%  before and after permutation to make sure that
%  that only the order and not the enties has changes;
%  exclude the sentinels, which are +-realmax and their
%  use may lead to rounding errors or overflows
% REVISION HISTORY:  
  len_a=length(a);
  list_2_index=1; % Start-Position for searching
  wrong_position_index_2=0;
  wrong_position_list_2=zeros(1,len_a);
% Initial sweep through the list
% write down the position where particles have to be interchanged
  while list_2_index<len_a-1
    if (a(list_2_index)>a(list_2_index+1))
      wrong_position_index_2=wrong_position_index_2+1;
      wrong_position_list_2(wrong_position_index_2)=list_2_index;
% If a particle position is marked for interchange, 
% Only the overnext particle has to be tested, because after
% interchange the new neighbors are checked of interchange in
% inc_sort
      list_2_index=list_2_index+2;
    else
      list_2_index=list_2_index+1;
    end
  end

%  search_depth=1;
%  num_perm=0;
  while wrong_position_index_2>0 % beginning of while-loop

%    num_perm=num_perm+wrong_position_index_2-1;
%    sum_a=sum(a(2:end-1));
    [wrong_position_index_1,wrong_position_list_1,a,b,c]=... 
    inc_sort_y2D(wrong_position_list_2,wrong_position_index_2,a,b,c);
%    if ~(abs(sum_a-sum(a(2:end-1)))*1e8<abs(sum_a))
%      error('Control sum of a_x-entries wrong 1')
%    end
%    search_depth=search_depth+1;
    if wrong_position_index_1==0 % verified ! 0, not 1
      break
    end

%    num_perm=num_perm+wrong_position_index_1-1;
%    sum_a=sum(a(2:end-1));
    [wrong_position_index_2,wrong_position_list_2,a,b,c]=... 
    inc_sort_y2D(wrong_position_list_1,wrong_position_index_1,a,b,c);
%    if ~(abs(sum_a-sum(a(2:end-1)))*1e8<abs(sum_a))
%      error('Control sum of a_y-entries wrong 2')
%    end
%    search_depth=search_depth+1;
    if wrong_position_index_1==0 % verified ! 0, not 1
      break
    end

  end   % end of while-loop
% end of permutations
return
