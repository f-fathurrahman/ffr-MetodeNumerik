disp('initializing capacitor updating coefficients');

for ind = 1:number_of_capacitors
    is = capacitors(ind).is;
    js = capacitors(ind).js;
    ks = capacitors(ind).ks;
    ie = capacitors(ind).ie;
    je = capacitors(ind).je;
    ke = capacitors(ind).ke;
    
    C = capacitors(ind).capacitance_per_component;
    
    switch (capacitors(ind).direction(1))
        case 'x'
            fi = create_linear_index_list(eps_r_x,is:ie-1,js:je,ks:ke);
            a_term = (2*C*dx)/(dy*dz);
            Cexe(fi) = ...
                (2*eps_0*eps_r_x(fi)-dt*sigma_e_x(fi)+a_term) ...
                    ./ (2*eps_0*eps_r_x(fi)+dt*sigma_e_x(fi)+a_term);
            Cexhz(fi)= (2*dt/dy)...
                ./ (2*eps_r_x(fi)*eps_0+dt*sigma_e_x(fi)+a_term);
            Cexhy(fi)= -(2*dt/dz) ...
                ./ (2*eps_r_x(fi)*eps_0 + dt*sigma_e_x(fi)+a_term);  
        case 'y'
            fi = create_linear_index_list(eps_r_y,is:ie,js:je-1,ks:ke);
            a_term = (2*C*dy)/(dz*dx);
            Ceye(fi) = ...
                (2*eps_0*eps_r_y(fi)-dt*sigma_e_y(fi)+a_term) ...
                    ./ (2*eps_0*eps_r_y(fi)+dt*sigma_e_y(fi)+a_term);
            Ceyhx(fi)= (2*dt/dz)...
                ./ (2*eps_r_y(fi)*eps_0+dt*sigma_e_y(fi)+a_term);
            Ceyhz(fi)= -(2*dt/dx) ...
                ./ (2*eps_r_y(fi)*eps_0 + dt*sigma_e_y(fi)+a_term);
        case 'z'
            fi = create_linear_index_list(eps_r_z,is:ie,js:je,ks:ke-1);
            a_term = (2*C*dz)/(dx*dy);
            Ceze(fi) = ...
                (2*eps_0*eps_r_z(fi)-dt*sigma_e_z(fi)+a_term) ...
                    ./ (2*eps_0*eps_r_z(fi)+dt*sigma_e_z(fi)+a_term);
            Cezhy(fi)= (2*dt/dx)...
                ./ (2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi)+a_term);
            Cezhx(fi)= -(2*dt/dy) ...
                ./ (2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi)+a_term);      
    end
    capacitors(ind).field_indices = fi;
end
