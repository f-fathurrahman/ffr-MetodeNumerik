function [node_indices] = get_node_indices(obj, fdtd_domain)
% convert coordinates to node indices on the FDTD grid

    if isfield(obj, 'min_x')
        node_indices.is = 1+round((obj.min_x - fdtd_domain.min_x)/fdtd_domain.dx); 
        node_indices.js = 1+round((obj.min_y - fdtd_domain.min_y)/fdtd_domain.dy); 
        node_indices.ks = 1+round((obj.min_z - fdtd_domain.min_z)/fdtd_domain.dz); 

        node_indices.ie = 1+round((obj.max_x - fdtd_domain.min_x)/fdtd_domain.dx); 
        node_indices.je = 1+round((obj.max_y - fdtd_domain.min_y)/fdtd_domain.dy); 
        node_indices.ke = 1+round((obj.max_z - fdtd_domain.min_z)/fdtd_domain.dz); 
    end

    if isfield(obj, 'x')
        node_indices.is = 1+round((obj.min_x - fdtd_domain.min_x)/fdtd_domain.dx); 
        node_indices.js = 1+round((obj.min_y - fdtd_domain.min_y)/fdtd_domain.dy); 
        node_indices.ks = 1+round((obj.min_z - fdtd_domain.min_z)/fdtd_domain.dz); 

        node_indices.ie = node_indices.is; 
        node_indices.je = node_indices.js; 
        node_indices.ke = node_indices.ks;         
    end