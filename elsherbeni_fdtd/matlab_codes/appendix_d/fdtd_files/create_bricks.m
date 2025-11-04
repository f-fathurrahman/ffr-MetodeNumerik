disp('creating bricks');

for ind = 1:number_of_bricks
    % convert brick end coordinates to node indices 
    ni = get_node_indices(bricks(ind), fdtd_domain);
    is = ni.is; js = ni.js; ks = ni.ks; 
    ie = ni.ie; je = ni.je; ke = ni.ke; 
    
    % assign material type of the brick to the cells
    material_3d_space (is:ie-1, js:je-1, ks:ke-1) ...
        = bricks(ind).material_type;
end
