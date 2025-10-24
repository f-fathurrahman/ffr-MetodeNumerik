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

