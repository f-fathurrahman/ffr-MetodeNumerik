function [node_indices] = get_node_indices(obj, fdtd_domain)
% convert coordinates to node indices on the FDTD grid

    if isfield(obj, 'min_x')
        [V, I] = min(abs(obj.min_x - fdtd_domain.node_coordinates_xe));
        node_indices.is = I(1); 
        [V, I] = min(abs(obj.min_y - fdtd_domain.node_coordinates_ye));
        node_indices.js = I(1); 
        [V, I] = min(abs(obj.min_z - fdtd_domain.node_coordinates_ze));
        node_indices.ks = I(1); 

        [V, I] = min(abs(obj.max_x - fdtd_domain.node_coordinates_xe));
        node_indices.ie = I(1); 
        [V, I] = min(abs(obj.max_y - fdtd_domain.node_coordinates_ye));
        node_indices.je = I(1); 
        [V, I] = min(abs(obj.max_z - fdtd_domain.node_coordinates_ze));
        node_indices.ke = I(1); 
    end

    if isfield(obj, 'x')
        [V, I] = min(abs(obj.x - fdtd_domain.node_coordinates_xe));
        node_indices.is = I(1); 
        [V, I] = min(abs(obj.y - fdtd_domain.node_coordinates_ye));
        node_indices.js = I(1); 
        [V, I] = min(abs(obj.z - fdtd_domain.node_coordinates_ze));
        node_indices.ks = I(1); 

        node_indices.ie = node_indices.is; 
        node_indices.je = node_indices.js; 
        node_indices.ke = node_indices.ks;         
    end
    
