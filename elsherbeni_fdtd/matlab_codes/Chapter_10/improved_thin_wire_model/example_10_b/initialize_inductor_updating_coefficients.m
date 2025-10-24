disp('initializing inductor updating coefficients');

for ind = 1:number_of_inductors
    is = inductors(ind).is;
    js = inductors(ind).js;
    ks = inductors(ind).ks;
    ie = inductors(ind).ie;
    je = inductors(ind).je;
    ke = inductors(ind).ke;

    L = inductors(ind).inductance_per_component;
        
    switch (inductors(ind).direction(1))
    case 'x'
        fi = create_linear_index_list(eps_r_x,is:ie-1,js:je,ks:ke);
        inductors(ind).Cexj = -(2*dt) ...
            ./ (2*eps_r_x(fi)*eps_0+dt*sigma_e_x(fi));      
        inductors(ind).Jix = zeros(size(fi));
        inductors(ind).Cjex = (dt*dx)/ (L*dy*dz);      
    case 'y'
        fi = create_linear_index_list(eps_r_y,is:ie,js:je-1,ks:ke);
        inductors(ind).Ceyj = -(2*dt) ...
            ./ (2*eps_r_y(fi)*eps_0+dt*sigma_e_y(fi));      
        inductors(ind).Jiy = zeros(size(fi));
        inductors(ind).Cjey = (dt*dy)/ (L*dz*dx);      
    case 'z'
        fi = create_linear_index_list(eps_r_z,is:ie,js:je,ks:ke-1);
        inductors(ind).Cezj = -(2*dt) ...
            ./ (2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi));      
        inductors(ind).Jiz = zeros(size(fi));
        inductors(ind).Cjez = (dt*dz)/ (L*dx*dy);      
    end
    inductors(ind).field_indices = fi;
end
