disp('initializing diode updating coefficients');

q = 1.602*1e-19; % charge of an electron
k = 1.38066e-23; % Boltzman constant, joule/kelvin 
T = 273+27;      % Kelvin; room temperature
I_0 = 1e-14;     % saturation current
    
for ind = 1:number_of_diodes
    is = diodes(ind).is;
    js = diodes(ind).js;
    ks = diodes(ind).ks;
    ie = diodes(ind).ie;
    je = diodes(ind).je;
    ke = diodes(ind).ke;
        
    if strcmp(diodes(ind).direction(2),'n')
        sgn = -1;
    else
        sgn = 1;
    end
    switch (diodes(ind).direction(1))
        case 'x'
            fi = create_linear_index_list(eps_r_x,is,js,ks);     
            diodes(ind).B = sgn*q*dx/(2*k*T);      
            diodes(ind).Cexd = ...
                -sgn*(2*dt*I_0/(dy*dz))*exp(diodes(ind).B) ...
                ./(2*eps_r_x(fi)*eps_0 + dt*sigma_e_x(fi));      
            diodes(ind).Exn = 0;      
        case 'y'
            fi = create_linear_index_list(eps_r_y,is,js,ks);     
            diodes(ind).B = sgn*q*dy/(2*k*T);      
            diodes(ind).Ceyd = ...
                -sgn*(2*dt*I_0/(dz*dx))*exp(diodes(ind).B) ...
                ./(2*eps_r_y(fi)*eps_0 + dt*sigma_e_y(fi));      
            diodes(ind).Eyn = 0;      
        case 'z'
            fi = create_linear_index_list(eps_r_z,is,js,ks);     
            diodes(ind).B = sgn*q*dz/(2*k*T);      
            diodes(ind).Cezd = ...
                -sgn*(2*dt*I_0/(dx*dy))*exp(diodes(ind).B) ...
                ./(2*eps_r_z(fi)*eps_0 + dt*sigma_e_z(fi));      
            diodes(ind).Ezn = 0;      
    end
    diodes(ind).field_indices = fi;
end
