disp('creating PEC plates on the material grid');

for ind = 1:number_of_bricks

    mtype     = bricks(ind).material_type;
    sigma_pec = material_types(mtype).sigma_e;

    % convert coordinates to node indices on the FDTD grid
    ni = get_node_indices(bricks(ind), fdtd_domain);
    is = ni.is; js = ni.js; ks = ni.ks; 
    ie = ni.ie; je = ni.je; ke = ni.ke; 
    
    % find the zero thickness bricks
    if (is == ie)
        sigma_e_y(is, js:je-1, ks:ke  ) = sigma_pec;
        sigma_e_z(is, js:je,   ks:ke-1) = sigma_pec;
    end
    if (js == je)
        sigma_e_z(is:ie,   js, ks:ke-1) = sigma_pec;
        sigma_e_x(is:ie-1, js, ks:ke  ) = sigma_pec;
    end
    if (ks == ke)
        sigma_e_x(is:ie-1, js:je,   ks) = sigma_pec;
        sigma_e_y(is:ie,   js:je-1, ks) = sigma_pec;
    end
end
