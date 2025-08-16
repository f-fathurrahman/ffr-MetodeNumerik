clear all
%   SYNTAX 
%   example_mesh_refinement
%   DESCRIPTION 
%   This script performs mesh refinement in a selected domain using
%   barycentric triangle subdivision. The input mesh must be 2 manifold
%
%   Low-Frequency Electromagnetic Modeling for Electrical and Biological
%   Systems Using MATLAB, Sergey N. Makarov, Gregory M. Noetscher, and Ara
%   Nazarian, Wiley, New York, 2015, 1st ed.

%%  Mesh (*.mat file)
FileName = 'VHP_Skull.mat';
load(FileName);
t = sort(t, 2);

%   Define the box where the mesh is refined
box = [-Inf Inf -Inf Inf 100 Inf];

%   Index into edges to be refined
edges = meshconnee(t);
AttachedTriangles = meshconnet(t, edges, 'manifold');

center  = 1/2*(P(edges(:, 1), :) + P(edges(:, 2), :));  
index1  = find(box(1)<center(:, 1)&center(:, 1)<box(2));
index2  = find(box(3)<center(:, 2)&center(:, 2)<box(4));
index3  = find(box(5)<center(:, 3)&center(:, 3)<box(6));
index   = intersect(index1, index2);
index   = intersect(index,  index3);

%   Nodes to be added up front
Nodes = center(index, :);
P = [Nodes; P];
t       = t + size(Nodes, 1);
edges   = edges + size(Nodes, 1);
%   Edges attached to every triangle
se      = meshconnte(t, edges);

%   Triangles to be added/removed
remove = [];
add    = [];
for m = 1:size(t, 1)
    temp1 = find(index==se(m, 1));
    temp2 = find(index==se(m, 2));
    temp3 = find(index==se(m, 3));
    node1 = intersect(edges(se(m, 1), :), edges(se(m, 3), :));
    node2 = intersect(edges(se(m, 1), :), edges(se(m, 2), :));
    node3 = intersect(edges(se(m, 2), :), edges(se(m, 3), :));
    if ~isempty(temp1)|~isempty(temp2)|~isempty(temp3)
        remove = [remove m];
        if temp1&temp2&temp3            
            add = [ [temp1 temp2 temp3];...
                    [temp1 temp2 node2];...
                    [temp1 temp3 node1];...               
                    [temp2 temp3 node3];...
                    add];
        end
        if temp1&temp2&(isempty(temp3))
            add = [ [temp1 temp2 node2];...
                    [temp1 temp2 node1];...               
                    [temp2 node1 node3];...
                    add];
        end
        if temp1&temp3&(isempty(temp2))
            add = [ [temp1 temp3 node1];...
                    [temp1 temp3 node3];...               
                    [temp1 node2 node3];...
                    add];
        end
        if temp2&temp3&(isempty(temp1))
            add = [ [temp2 temp3 node3];...
                    [temp2 temp3 node1];...               
                    [temp2 node1 node2];...
                    add];
        end
        if temp1&(isempty(temp2))&(isempty(temp3))
            add = [ [temp1 node1 node3];...    
                    [temp1 node2 node3];...
                    add];
        end
        if temp2&(isempty(temp1))&(isempty(temp3))
            add = [ [temp2 node1 node2];...    
                    [temp2 node1 node3];...
                    add];
        end
        if temp3&(isempty(temp1))&(isempty(temp2))
            add = [ [temp3 node1 node2];...    
                    [temp3 node2 node3];...
                    add];   
        end
    end
end

t(remove, :)    = [];
t               = [t; add];
t               = sort(t, 2);

X = reshape(P(t', 1),[3, size(t, 1)]);
Y = reshape(P(t', 2),[3, size(t, 1)]);
Z = reshape(P(t', 3),[3, size(t, 1)]);
patch(X, Y, Z, [1 0.75 0.65], 'EdgeColor', 'k', 'FaceAlpha', 1.0, 'AmbientStrength', 0.75);

light('Position', [0 1 0], 'Style', 'local'); axis 'equal';  axis 'tight'; view(0, 0); grid on;
xlabel('x, mm'); ylabel('y, mm'); zlabel('z, mm');

save(strcat(FileName(1:end-4),'_mod'), 'P', 't');







