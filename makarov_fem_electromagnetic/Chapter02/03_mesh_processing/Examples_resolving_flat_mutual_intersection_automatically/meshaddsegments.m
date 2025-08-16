function si13 = meshaddsegments(st1, st2, edges_t1, edges_t2, si11, si12, si21, si22, SI11, SI21);
%   SYNTAX 
%   si13 = meshaddsegments(st1, st2, edges_t1, edges_t2, si11, si12, si21, si22, SI11, SI21);
%   DESCRIPTION 
%   For every triangle of intersected mesh #1, this function finds extra segments (node pairs) to be added
%   Duplicated points will be present
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

    si13 = cell(st1, 1);            %   extra segments
    for m = SI11                    %   loop over suspicious triangles of mesh #1
        temp = [];                  %   collects all intersecting edges of the slave mesh(sorted) for triangle m of the master mesh       
        for n = 1:size(si11, 1)
            if sum(si11{n}==m)>0    %   edge n of mesh #2 intersects triangle m of mesh #1 (and maybe other triangles) 
                temp  = [temp n];
            end
        end 
        if ~isempty(temp)       
            for n = SI21            %   loop over suspicious triangles of mesh #2
                temp2 = intersect(temp, edges_t2(n ,:)); 
                %   Case #1: exactly two edges in temp2 of the slave triangle n intersect the
                %   master triangle m. One intersection edge has to be added. 
                if length(temp2)==2 % full triangle found
                    index1 = find(si11{temp2(1)}==m);
                    point1 = si12{temp2(1)}(index1, :);
                    index2 = find(si11{temp2(2)}==m);
                    point2 = si12{temp2(2)}(index2, :);                
                    si13{m} = [si13{m}; point1; point2];
                end
                %   Case #2: only one edge in temp2 of the slave triangle n intersects the 
                %   master triangle m. Then, an edge of the master triangle must
                %   intersect this slave triangle. One intersection edge has to be added.             
                if length(temp2)==1 % "half" triangle found
                    index1 = find(si11{temp2}==m);
                    point1 = si12{temp2}(index1, :);
                    for p = 1:3
                        edgemaster = edges_t1(m, p);
                        index2 = find(si21{edgemaster}==n);
                        if ~isempty(index2)
                            point2 = si22{edgemaster}(index2, :);
                            si13{m} = [si13{m}; point1; point2];
                        end
                    end                                
                end
            end
        end 
        %   Case #3: two edges of the master triangle intersect the same slave
        %   triangle. One intersection edge has to be added. 
        trimaster1 = si21{edges_t1(m, 1)}';     %   row of intersected slave triangles (numbers)
        trimaster2 = si21{edges_t1(m, 2)}';     %   row of intersected slave triangles (numbers)
        trimaster3 = si21{edges_t1(m, 3)}';     %   row of intersected slave triangles (numbers)
        index1 = intersect(trimaster1, trimaster2); % numbers of intersected slave triangles
        if ~isempty(index1)
            for p = 1:length(index1)
                point1 = si22{edges_t1(m, 1)}(find(trimaster1==index1(p)), :);
                point2 = si22{edges_t1(m, 2)}(find(trimaster2==index1(p)), :);
                si13{m} = [si13{m}; point1; point2];
            end
        end
        index1 = intersect(trimaster1, trimaster3); % numbers of intersected slave triangles
        if ~isempty(index1)
            for p = 1:length(index1)
                point1 = si22{edges_t1(m, 1)}(find(trimaster1==index1(p)), :);
                point2 = si22{edges_t1(m, 3)}(find(trimaster3==index1(p)), :);
                si13{m} = [si13{m}; point1; point2];
            end
        end        
        index1 = intersect(trimaster2, trimaster3); % numbers of intersected slave triangles
        if ~isempty(index1)
            for p = 1:length(index1)
                point1 = si22{edges_t1(m, 2)}(find(trimaster2==index1(p)), :);
                point2 = si22{edges_t1(m, 3)}(find(trimaster3==index1(p)), :);
                si13{m} = [si13{m}; point1; point2];
            end
        end
    end
end