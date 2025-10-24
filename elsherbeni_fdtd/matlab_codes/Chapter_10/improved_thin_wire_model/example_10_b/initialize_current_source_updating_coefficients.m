disp('initializing current source updating coefficients');

for ind = 1:number_of_current_sources
    is = current_sources(ind).is;
    js = current_sources(ind).js;
    ks = current_sources(ind).ks;
    ie = current_sources(ind).ie;
    je = current_sources(ind).je;
    ke = current_sources(ind).ke;
    
    R = current_sources(ind).resistance_per_component;
    
    switch (current_sources(ind).direction(1))
    case 'x'
        fi = create_linear_index_list(eps_r_x,is:ie-1,js:je,ks:ke);
        a_term = (dt*dx)/(R*dy*dz);
        Cexe(fi) = ...
            (2*eps_0*eps_r_x(fi)-dt*sigma_e_x(fi)-a_term) ...
                ./ (2*eps_0*eps_r_x(fi)+dt*sigma_e_x(fi)+a_term);
        Cexhz(fi)= (2*dt/dy)...
            ./ (2*eps_r_x(fi)*eps_0+dt*sigma_e_x(fi)+a_term);
        Cexhy(fi)= -(2*dt/dz) ...
            ./ (2*eps_r_x(fi)*eps_0 + dt*sigma_e_x(fi)+a_term);  
        current_sources(ind).Cexs = -(2*dt/(dy*dz)) ...
            ./(2*eps_r_x(fi)*eps_0+dt*sigma_e_x(fi)+a_term);
    case 'y'
        fi = create_linear_index_list(eps_r_y,is:ie,js:je-1,ks:ke);
        a_term = (dt*dy)/(R*dz*dx);
        Ceye(fi) = ...
            (2*eps_0*eps_r_y(fi)-dt*sigma_e_y(fi)-a_term) ...
                ./ (2*eps_0*eps_r_y(fi)+dt*sigma_e_y(fi)+a_term);
        Ceyhx(fi)= (2*dt/dz)...
            ./ (2*eps_r_y(fi)*eps_0+dt*sigma_e_y(fi)+a_term);
        Ceyhz(fi)= -(2*dt/dx) ...
            ./ (2*eps_r_y(fi)*eps_0 + dt*sigma_e_y(fi)+a_term);
        current_sources(ind).Ceys = -(2*dt/(dz*dx)) ...
            ./(2*eps_r_y(fi)*eps_0+dt*sigma_e_y(fi)+a_term);
    case 'z'
        fi = create_linear_index_list(eps_r_z,is:ie,js:je,ks:ke-1);
        a_term = (dt*dz)/(R*dx*dy);
        Ceze(fi) = ...
            (2*eps_0*eps_r_z(fi)-dt*sigma_e_z(fi)-a_term) ...
                ./ (2*eps_0*eps_r_z(fi)+dt*sigma_e_z(fi)+a_term);
        Cezhy(fi)= (2*dt/dx)...
            ./ (2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi)+a_term);
        Cezhx(fi)= -(2*dt/dy) ...
            ./ (2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi)+a_term);      
        current_sources(ind).Cezs = -(2*dt/(dx*dy)) ...
            ./(2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi)+a_term);
    end
    current_sources(ind).field_indices = fi;
end
