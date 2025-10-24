disp('initializing sources and lumped element components');

number_of_voltage_sources   = size(voltage_sources,2);
number_of_current_sources   = size(current_sources,2);
number_of_resistors   = size(resistors,2);
number_of_inductors   = size(inductors,2);
number_of_capacitors  = size(capacitors,2);
number_of_diodes      = size(diodes,2);

% initialize waveforms
initialize_waveforms;

% voltage sources
for ind = 1:number_of_voltage_sources
    is = round((voltage_sources(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((voltage_sources(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((voltage_sources(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((voltage_sources(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((voltage_sources(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((voltage_sources(ind).max_z - fdtd_domain.min_z)/dz)+1;
    voltage_sources(ind).is = is;
    voltage_sources(ind).js = js;
    voltage_sources(ind).ks = ks;
    voltage_sources(ind).ie = ie;
    voltage_sources(ind).je = je;
    voltage_sources(ind).ke = ke;

    switch (voltage_sources(ind).direction(1))
    case 'x'
        n_fields = ie - is;
        r_magnitude_factor = (1 + je - js) * (1 + ke - ks) / (ie - is); 
    case 'y'
        n_fields = je - js;
        r_magnitude_factor = (1 + ie - is) * (1 + ke - ks) / (je - js); 
    case 'z'
        n_fields = ke - ks;
        r_magnitude_factor = (1 + ie - is) * (1 + je - js) / (ke - ks); 
    end
    if strcmp(voltage_sources(ind).direction(2),'n')
        v_magnitude_factor =  ...
            -1*voltage_sources(ind).magnitude/n_fields;
    else
        v_magnitude_factor =  ...
            1*voltage_sources(ind).magnitude/n_fields;
    end
    voltage_sources(ind).resistance_per_component = ...
        r_magnitude_factor * voltage_sources(ind).resistance;
    
    % copy waveform of the waveform type to waveform of the source 
    wt_str = voltage_sources(ind).waveform_type;
    wi_str = num2str(voltage_sources(ind).waveform_index);
    eval_str = ['a_waveform = waveforms.' ...
        wt_str '(' wi_str ').waveform;'];
    eval(eval_str);
    voltage_sources(ind).voltage_per_e_field = ...
        v_magnitude_factor * a_waveform;
    voltage_sources(ind).waveform = ...
        v_magnitude_factor * a_waveform * n_fields;
end

% current sources
for ind = 1:number_of_current_sources
    is = round((current_sources(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((current_sources(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((current_sources(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((current_sources(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((current_sources(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((current_sources(ind).max_z - fdtd_domain.min_z)/dz)+1;
    current_sources(ind).is = is;
    current_sources(ind).js = js;
    current_sources(ind).ks = ks;
    current_sources(ind).ie = ie;
    current_sources(ind).je = je;
    current_sources(ind).ke = ke;

    switch (current_sources(ind).direction(1))
    case 'x'
        n_fields = (1 + je - js) * (1 + ke - ks);
        r_magnitude_factor = (1 + je - js) * (1 + ke - ks) / (ie - is); 
    case 'y'
        n_fields = (1 + ie - is) * (1 + ke - ks);
        r_magnitude_factor = (1 + ie - is) * (1 + ke - ks) / (je - js); 
    case 'z'
        n_fields = (1 + ie - is) * (1 + je - js);
        r_magnitude_factor = (1 + ie - is) * (1 + je - js) / (ke - ks); 
    end
    if strcmp(current_sources(ind).direction(2),'n')
        i_magnitude_factor = ...
            -1*current_sources(ind).magnitude/n_fields;
    else
        i_magnitude_factor =  ...
            1*current_sources(ind).magnitude/n_fields;
    end
    current_sources(ind).resistance_per_component = ...
        r_magnitude_factor * current_sources(ind).resistance;

    % copy waveform of the waveform type to waveform of the source 
    wt_str = current_sources(ind).waveform_type;
    wi_str = num2str(current_sources(ind).waveform_index);
    eval_str = ['a_waveform = waveforms.' ...
        wt_str '(' wi_str ').waveform;'];
    eval(eval_str);
    current_sources(ind).current_per_e_field = ...
        i_magnitude_factor * a_waveform;
    current_sources(ind).waveform = ...
        i_magnitude_factor * a_waveform * n_fields;
end

% resistors
for ind = 1:number_of_resistors
    is = round((resistors(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((resistors(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((resistors(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((resistors(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((resistors(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((resistors(ind).max_z - fdtd_domain.min_z)/dz)+1;
    resistors(ind).is = is;
    resistors(ind).js = js;
    resistors(ind).ks = ks;
    resistors(ind).ie = ie;
    resistors(ind).je = je;
    resistors(ind).ke = ke;
    switch (resistors(ind).direction)
    case 'x'
        r_magnitude_factor = (1 + je - js) * (1 + ke - ks) / (ie - is); 
    case 'y'
        r_magnitude_factor = (1 + ie - is) * (1 + ke - ks) / (je - js); 
    case 'z'
        r_magnitude_factor = (1 + ie - is) * (1 + je - js) / (ke - ks); 
    end
    resistors(ind).resistance_per_component = ...
        r_magnitude_factor * resistors(ind).resistance;
end

% inductors
for ind = 1:number_of_inductors
    is = round((inductors(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((inductors(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((inductors(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((inductors(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((inductors(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((inductors(ind).max_z - fdtd_domain.min_z)/dz)+1;
    inductors(ind).is = is;
    inductors(ind).js = js;
    inductors(ind).ks = ks;
    inductors(ind).ie = ie;
    inductors(ind).je = je;
    inductors(ind).ke = ke;
    switch (inductors(ind).direction)
    case 'x'
        l_magnitude_factor = (1 + je - js) * (1 + ke - ks) / (ie - is); 
    case 'y'
        l_magnitude_factor = (1 + ie - is) * (1 + ke - ks) / (je - js); 
    case 'z'
        l_magnitude_factor = (1 + ie - is) * (1 + je - js) / (ke - ks); 
    end
    inductors(ind).inductance_per_component = ...
        l_magnitude_factor * inductors(ind).inductance;
end

% capacitors
for ind = 1:number_of_capacitors
    is = round((capacitors(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((capacitors(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((capacitors(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((capacitors(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((capacitors(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((capacitors(ind).max_z - fdtd_domain.min_z)/dz)+1;
    capacitors(ind).is = is;
    capacitors(ind).js = js;
    capacitors(ind).ks = ks;
    capacitors(ind).ie = ie;
    capacitors(ind).je = je;
    capacitors(ind).ke = ke;
    switch (capacitors(ind).direction)
    case 'x'
        c_magnitude_factor = (ie - is) ...
            / ((1 + je - js) * (1 + ke - ks)); 
    case 'y'
        c_magnitude_factor = (je - js) ...
            / ((1 + ie - is) * (1 + ke - ks)); 
    case 'z'
        c_magnitude_factor = (ke - ks) ...
            / ((1 + ie - is) * (1 + je - js)); 
    end
    capacitors(ind).capacitance_per_component = ...
        c_magnitude_factor * capacitors(ind).capacitance;
end

sigma_pec = material_types(material_type_index_pec).sigma_e;
 
% diodes
for ind = 1:number_of_diodes
    is = round((diodes(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((diodes(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((diodes(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((diodes(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((diodes(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((diodes(ind).max_z - fdtd_domain.min_z)/dz)+1;
    diodes(ind).is = is;
    diodes(ind).js = js;
    diodes(ind).ks = ks;
    diodes(ind).ie = ie;
    diodes(ind).je = je;
    diodes(ind).ke = ke;

    switch (diodes(ind).direction(1))
        case 'x'
            sigma_e_x(is+1:ie-1,js,ks) = sigma_pec;
        case 'y'
            sigma_e_y(is,js+1:je-1,ks) = sigma_pec;
        case 'z'
            sigma_e_z(is,js,ks+1:ke-1) = sigma_pec;
    end
end
