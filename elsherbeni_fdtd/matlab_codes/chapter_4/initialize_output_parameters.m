disp('initializing the output parameters');

number_of_sampled_electric_fields  = size(sampled_electric_fields,2);
number_of_sampled_magnetic_fields  = size(sampled_magnetic_fields,2);
number_of_sampled_voltages  = size(sampled_voltages,2);
number_of_sampled_currents  = size(sampled_currents,2);

% initialize sampled electric field terms
for ind=1:number_of_sampled_electric_fields  
    is = round((sampled_electric_fields(ind).x ...
        - fdtd_domain.min_x)/dx)+1;
    js = round((sampled_electric_fields(ind).y ...
        - fdtd_domain.min_y)/dy)+1;
    ks = round((sampled_electric_fields(ind).z ...
        - fdtd_domain.min_z)/dz)+1;
    sampled_electric_fields(ind).is = is;
    sampled_electric_fields(ind).js = js;
    sampled_electric_fields(ind).ks = ks;
    sampled_electric_fields(ind).sampled_value = ...
        zeros(1, number_of_time_steps);
    sampled_electric_fields(ind).time = ...
        ([1:number_of_time_steps])*dt;
end

% initialize sampled magnetic field terms
for ind=1:number_of_sampled_magnetic_fields  
    is = round((sampled_magnetic_fields(ind).x ...
        - fdtd_domain.min_x)/dx)+1;
    js = round((sampled_magnetic_fields(ind).y ...
        - fdtd_domain.min_y)/dy)+1;
    ks = round((sampled_magnetic_fields(ind).z ...
        - fdtd_domain.min_z)/dz)+1;
    sampled_magnetic_fields(ind).is = is;
    sampled_magnetic_fields(ind).js = js;
    sampled_magnetic_fields(ind).ks = ks;
    sampled_magnetic_fields(ind).sampled_value = ...
        zeros(1, number_of_time_steps);
    sampled_magnetic_fields(ind).time = ...
        ([1:number_of_time_steps]-0.5)*dt;
end

% initialize sampled voltage terms
for ind=1:number_of_sampled_voltages  
    is = round((sampled_voltages(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((sampled_voltages(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((sampled_voltages(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((sampled_voltages(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((sampled_voltages(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((sampled_voltages(ind).max_z - fdtd_domain.min_z)/dz)+1;
    sampled_voltages(ind).is = is;
    sampled_voltages(ind).js = js;
    sampled_voltages(ind).ks = ks;
    sampled_voltages(ind).ie = ie;
    sampled_voltages(ind).je = je;
    sampled_voltages(ind).ke = ke;
    sampled_voltages(ind).sampled_value = ...
                        zeros(1, number_of_time_steps);

    switch (sampled_voltages(ind).direction(1))
    case 'x'
        fi = create_linear_index_list(Ex,is:ie-1,js:je,ks:ke);
        sampled_voltages(ind).Csvf = -dx/((je-js+1)*(ke-ks+1));
    case 'y'
        fi = create_linear_index_list(Ey,is:ie,js:je-1,ks:ke);
        sampled_voltages(ind).Csvf = -dy/((ke-ks+1)*(ie-is+1));
    case 'z'
        fi = create_linear_index_list(Ez,is:ie,js:je,ks:ke-1);
        sampled_voltages(ind).Csvf = -dz/((ie-is+1)*(je-js+1));
    end    
    if strcmp(sampled_voltages(ind).direction(2),'n')
        sampled_voltages(ind).Csvf =  ...
            -1 * sampled_voltages(ind).Csvf;
    end
    sampled_voltages(ind).field_indices = fi;
    sampled_voltages(ind).time = ([1:number_of_time_steps])*dt;
end

% initialize sampled current terms
for ind=1:number_of_sampled_currents  
    is = round((sampled_currents(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((sampled_currents(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((sampled_currents(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((sampled_currents(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((sampled_currents(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((sampled_currents(ind).max_z - fdtd_domain.min_z)/dz)+1;
    sampled_currents(ind).is = is;
    sampled_currents(ind).js = js;
    sampled_currents(ind).ks = ks;
    sampled_currents(ind).ie = ie;
    sampled_currents(ind).je = je;
    sampled_currents(ind).ke = ke;
    sampled_currents(ind).sampled_value = ...
                        zeros(1, number_of_time_steps);
    sampled_currents(ind).time =([1:number_of_time_steps]-0.5)*dt;
end
