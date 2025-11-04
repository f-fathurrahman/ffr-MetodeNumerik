disp('initializing diode updating coefficients');

q = 1.602*1e-19; % charge of an electron
k = 1.38066e-23; % Boltzman constant, joule/kelvin 
T = 273+27;      % Kelvin; room temperature
I_0 = 1e-14;     % saturation current
    
for ind = 1:number_of_diodes
    ni = diodes(ind).node_indices;
    is = ni.is; js = ni.js; ks = ni.ks; 
    ie = ni.ie; je = ni.je; ke = ni.ke; 
        
    if strcmp(diodes(ind).direction(2),'n')
        sgn = -1;
    else
        sgn = 1;
    end
    switch (diodes(ind).direction(1))
        case 'x'
            fi = create_linear_index_list(eps_r_x,is,js,ks);

            DX = cell_sizes_xe(is);
            DY = cell_sizes_yh(js);
            DZ = cell_sizes_zh(ks);
            
            diodes(ind).B = sgn*q*DX/(2*k*T);      
            diodes(ind).Cexd = -sgn*(2*dt*I_0/(DY*DZ)) ...
                ./(2*eps_r_x(fi)*eps_0 + dt*sigma_e_x(fi));      
            diodes(ind).Exn = 0;      
        case 'y'
            fi = create_linear_index_list(eps_r_y,is,js,ks);     

            DX = cell_sizes_xh(is);
            DY = cell_sizes_ye(js);
            DZ = cell_sizes_zh(ks);
            
            diodes(ind).B = sgn*q*DY/(2*k*T);      
            diodes(ind).Ceyd = -sgn*(2*dt*I_0/(DZ*DX)) ...
                ./(2*eps_r_y(fi)*eps_0 + dt*sigma_e_y(fi));      
            diodes(ind).Eyn = 0;      
        case 'z'
            fi = create_linear_index_list(eps_r_z,is,js,ks);     

            DX = cell_sizes_xh(is);
            DY = cell_sizes_yh(js);
            DZ = cell_sizes_ze(ks);
            
            diodes(ind).B = sgn*q*DZ/(2*k*T);      
            diodes(ind).Cezd = -sgn*(2*dt*I_0/(DX*DY)) ...
                ./(2*eps_r_z(fi)*eps_0 + dt*sigma_e_z(fi));      
            diodes(ind).Ezn = 0;      
    end
    diodes(ind).field_indices = fi;
end
