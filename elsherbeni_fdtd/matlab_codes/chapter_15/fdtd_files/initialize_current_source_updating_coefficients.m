disp('initializing current source updating coefficients');

for ind = 1:number_of_current_sources
    ni = current_sources(ind).node_indices;
    is = ni.is; js = ni.js; ks = ni.ks; 
    ie = ni.ie; je = ni.je; ke = ni.ke; 
    
    R = current_sources(ind).resistance_per_component;
    
    switch (current_sources(ind).direction(1))
    case 'x'
        fi = create_linear_index_list(eps_r_x,is:ie-1,js:je,ks:ke);
        
        csx = cell_sizes_xe(is:ie-1);
        csy = cell_sizes_yh(js:je);
        csz = cell_sizes_zh(ks:ke);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DX = reshape(DX,[],1);
        DY = reshape(DY,[],1);
        DZ = reshape(DZ,[],1);

        a_term = (dt*DX)./(R*DY.*DZ);
        Cexe(fi) = ...
            (2*eps_0*eps_r_x(fi)-dt*sigma_e_x(fi)-a_term) ...
                ./ (2*eps_0*eps_r_x(fi)+dt*sigma_e_x(fi)+a_term);
        Cexhz(fi)= (2*dt./DY)...
            ./ (2*eps_r_x(fi)*eps_0+dt*sigma_e_x(fi)+a_term);
        Cexhy(fi)= -(2*dt./DZ) ...
            ./ (2*eps_r_x(fi)*eps_0 + dt*sigma_e_x(fi)+a_term);  
        current_sources(ind).Cexs = -(2*dt./(DY.*DZ)) ...
            ./(2*eps_r_x(fi)*eps_0+dt*sigma_e_x(fi)+a_term);
    case 'y'
        fi = create_linear_index_list(eps_r_y,is:ie,js:je-1,ks:ke);
        
        csx = cell_sizes_xh(is:ie);
        csy = cell_sizes_ye(js:je-1);
        csz = cell_sizes_zh(ks:ke);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DX = reshape(DX,[],1);
        DY = reshape(DY,[],1);
        DZ = reshape(DZ,[],1);

        a_term = (dt*DY)./(R*DZ.*DX);
        Ceye(fi) = ...
            (2*eps_0*eps_r_y(fi)-dt*sigma_e_y(fi)-a_term) ...
                ./ (2*eps_0*eps_r_y(fi)+dt*sigma_e_y(fi)+a_term);
        Ceyhx(fi)= (2*dt./DZ)...
            ./ (2*eps_r_y(fi)*eps_0+dt*sigma_e_y(fi)+a_term);
        Ceyhz(fi)= -(2*dt./DX) ...
            ./ (2*eps_r_y(fi)*eps_0 + dt*sigma_e_y(fi)+a_term);
        current_sources(ind).Ceys = -(2*dt./(DZ.*DX)) ...
            ./(2*eps_r_y(fi)*eps_0+dt*sigma_e_y(fi)+a_term);
    case 'z'
        fi = create_linear_index_list(eps_r_z,is:ie,js:je,ks:ke-1);

        csx = cell_sizes_xh(is:ie);
        csy = cell_sizes_yh(js:je);
        csz = cell_sizes_ze(ks:ke-1);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DX = reshape(DX,[],1);
        DY = reshape(DY,[],1);
        DZ = reshape(DZ,[],1);

        a_term = (dt*DZ)./(R*DX.*DY);
        Ceze(fi) = ...
            (2*eps_0*eps_r_z(fi)-dt*sigma_e_z(fi)-a_term) ...
                ./ (2*eps_0*eps_r_z(fi)+dt*sigma_e_z(fi)+a_term);
        Cezhy(fi)= (2*dt./DX)...
            ./ (2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi)+a_term);
        Cezhx(fi)= -(2*dt./DY) ...
            ./ (2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi)+a_term);      
        current_sources(ind).Cezs = -(2*dt./(DX.*DY)) ...
            ./(2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi)+a_term);
    end
    current_sources(ind).field_indices = fi;
end
