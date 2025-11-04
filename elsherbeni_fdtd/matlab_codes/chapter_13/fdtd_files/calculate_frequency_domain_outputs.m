disp('generating frequency domain outputs');

frequency_array = frequency_domain.frequencies;

% sampled electric fields in frequency domain
for ind=1:number_of_sampled_electric_fields
    x = sampled_electric_fields(ind).sampled_value;
    time_shift = 0;
    [X] = time_to_frequency_domain(x, dt, frequency_array, time_shift);
    sampled_electric_fields(ind).frequency_domain_value = X;
    sampled_electric_fields(ind).frequencies = frequency_array;
end

% sampled magnetic fields in frequency domain
for ind=1:number_of_sampled_magnetic_fields
    x = sampled_magnetic_fields(ind).sampled_value;
    time_shift = -dt/2;
    [X] = time_to_frequency_domain(x, dt, frequency_array, time_shift);
    sampled_magnetic_fields(ind).frequency_domain_value = X;
    sampled_magnetic_fields(ind).frequencies = frequency_array;
end

% sampled voltages in frequency domain
for ind=1:number_of_sampled_voltages
    x = sampled_voltages(ind).sampled_value;
    time_shift = 0;
    [X] = time_to_frequency_domain(x, dt, frequency_array, time_shift);
    sampled_voltages(ind).frequency_domain_value = X;
    sampled_voltages(ind).frequencies = frequency_array;
end

% sampled currents in frequency domain
for ind=1:number_of_sampled_currents
    x = sampled_currents(ind).sampled_value;
    time_shift = -dt/2;
    [X] = time_to_frequency_domain(x, dt, frequency_array, time_shift);
    sampled_currents(ind).frequency_domain_value = X;
    sampled_currents(ind).frequencies = frequency_array;
end

% voltage sources in frequency domain
for ind=1:number_of_voltage_sources
    x = voltage_sources(ind).waveform;
    time_shift = 0;
    [X] = time_to_frequency_domain(x, dt, frequency_array, time_shift);
    voltage_sources(ind).frequency_domain_value = X;
    voltage_sources(ind).frequencies = frequency_array;
end

% current sources in frequency domain
for ind=1:number_of_current_sources
    x = current_sources(ind).waveform;
    time_shift = 0;
    [X] = time_to_frequency_domain(x, dt, frequency_array, time_shift);
    current_sources(ind).frequency_domain_value = X;
    current_sources(ind).frequencies = frequency_array;
end

% calculation of S-parameters
% calculate incident and reflected power waves
for ind=1:number_of_ports
   svi = ports(ind).sampled_voltage_index; 
   sci = ports(ind).sampled_current_index; 
   Z = ports(ind).impedance;
   V = sampled_voltages(svi).frequency_domain_value;
   I = sampled_currents(sci).frequency_domain_value;
   ports(ind).a = 0.5*(V+Z.*I)./sqrt(real(Z));
   ports(ind).b = 0.5*(V-conj(Z).*I)./sqrt(real(Z));   
   ports(ind).frequencies = frequency_array;
end

% calculate the S-parameters
for ind=1:number_of_ports
    if ports(ind).is_source_port == true
        for oind=1:number_of_ports
            ports(ind).S(oind).values = ports(oind).b ./ ports(ind).a;
        end
    end
end
