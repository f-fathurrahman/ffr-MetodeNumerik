disp('initializing sources and lumped element components');

number_of_voltage_sources   = 0;
number_of_current_sources   = 0;
number_of_resistors   = 0;
number_of_inductors   = 0;
number_of_capacitors  = 0;
number_of_diodes      = 0;

if exist('voltage_sources','var')
    number_of_voltage_sources = size(voltage_sources,2);
end
if exist('current_sources','var')
    number_of_current_sources = size(current_sources,2);
end
if exist('resistors','var')
    number_of_resistors = size(resistors,2);
end
if exist('inductors','var')
    number_of_inductors = size(inductors,2);
end
if exist('capacitors','var')
    number_of_capacitors = size(capacitors,2);
end
if exist('diodes','var')
    number_of_diodes = size(diodes,2);
end

% voltage sources
for ind = 1:number_of_voltage_sources
    node_indices = get_node_indices(voltage_sources(ind), fdtd_domain);
    is = node_indices.is; js = node_indices.js; ks = node_indices.ks; 
    ie = node_indices.ie; je = node_indices.je; ke = node_indices.ke; 
    
    voltage_sources(ind).node_indices = node_indices;

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
        v_magnitude_factor =  -1*voltage_sources(ind).magnitude/n_fields;
    else
        v_magnitude_factor =  1*voltage_sources(ind).magnitude/n_fields;
    end
    voltage_sources(ind).resistance_per_component = ...
        r_magnitude_factor * voltage_sources(ind).resistance;
    
    current_object = voltage_sources(ind);
    calculate_waveform;
    voltage_sources(ind).waveform = current_object.waveform;
    
    voltage_sources(ind).voltage_per_e_field = ...
        v_magnitude_factor * current_object.waveform;
end

% current sources
for ind = 1:number_of_current_sources
    node_indices = get_node_indices(current_sources(ind), fdtd_domain);
    is = node_indices.is; js = node_indices.js; ks = node_indices.ks; 
    ie = node_indices.ie; je = node_indices.je; ke = node_indices.ke; 
    
    current_sources(ind).node_indices = node_indices;

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
        i_magnitude_factor = -1*current_sources(ind).magnitude/n_fields;
    else
        i_magnitude_factor =  current_sources(ind).magnitude/n_fields;
    end
    current_sources(ind).resistance_per_component = ...
        r_magnitude_factor * current_sources(ind).resistance;

    current_object = current_sources(ind);
    calculate_waveform;
    current_sources(ind).waveform = current_object.waveform;
    current_sources(ind).current_per_e_field = ...
        i_magnitude_factor * current_object.waveform;
end

% resistors
for ind = 1:number_of_resistors

    node_indices = get_node_indices(resistors(ind), fdtd_domain);
    is = node_indices.is; js = node_indices.js; ks = node_indices.ks; 
    ie = node_indices.ie; je = node_indices.je; ke = node_indices.ke; 
    
    resistors(ind).node_indices = node_indices;
    
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
    
    node_indices = get_node_indices(inductors(ind), fdtd_domain);
    is = node_indices.is; js = node_indices.js; ks = node_indices.ks; 
    ie = node_indices.ie; je = node_indices.je; ke = node_indices.ke; 
    
    inductors(ind).node_indices = node_indices;
    
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
    
    node_indices = get_node_indices(capacitors(ind), fdtd_domain);
    is = node_indices.is; js = node_indices.js; ks = node_indices.ks; 
    ie = node_indices.ie; je = node_indices.je; ke = node_indices.ke; 
    
    capacitors(ind).node_indices = node_indices;

    switch (capacitors(ind).direction)
        case 'x'
            c_magnitude_factor = (ie - is) / ((1 + je - js) * (1 + ke - ks)); 
        case 'y'
            c_magnitude_factor = (je - js) / ((1 + ie - is) * (1 + ke - ks)); 
        case 'z'
            c_magnitude_factor = (ke - ks) / ((1 + ie - is) * (1 + je - js)); 
    end
    capacitors(ind).capacitance_per_component = ...
        c_magnitude_factor * capacitors(ind).capacitance;
end
 
% diodes
for ind = 1:number_of_diodes
    
    node_indices = get_node_indices(diodes(ind), fdtd_domain);
    is = node_indices.is; js = node_indices.js; ks = node_indices.ks; 
    ie = node_indices.ie; je = node_indices.je; ke = node_indices.ke; 
    
    diodes(ind).node_indices = node_indices;

    switch (diodes(ind).direction(1))
        case 'x'
            sigma_e_x(is+1:ie-1,js,ks) = sigma_pec;
        case 'y'
            sigma_e_y(is,js+1:je-1,ks) = sigma_pec;
        case 'z'
            sigma_e_z(is,js,ks+1:ke-1) = sigma_pec;
    end
end

% initialize incident plane wave
if exist('incident_plane_wave','var')
    if isfield(incident_plane_wave,'enable')
        if incident_plane_wave.enable 
            initialize_incident_plane_wave;            
        end
    else
        incident_plane_wave.enable = true;
        initialize_incident_plane_wave;
    end
else
    incident_plane_wave.enable = false;
end
