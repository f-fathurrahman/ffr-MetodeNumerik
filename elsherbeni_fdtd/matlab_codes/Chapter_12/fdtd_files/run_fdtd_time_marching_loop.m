disp (['Starting the time marching loop']);
disp(['Total number of time steps : ' ...
    num2str(number_of_time_steps)]); 

start_time = cputime; 
current_time = 0;

for time_step = 1:number_of_time_steps  
    update_incident_fields;
    update_magnetic_fields;
    update_magnetic_field_CPML_ABC;
    capture_sampled_magnetic_fields;
    capture_sampled_currents;
    update_electric_fields;
    update_electric_field_CPML_ABC;
    update_voltage_sources; 
    update_current_sources; 
    update_inductors; 
    update_diodes; 
    capture_sampled_electric_fields;
    capture_sampled_voltages;
    calculate_JandM; 
    display_sampled_parameters;
end                                 

end_time = cputime;
total_time_in_minutes = (end_time - start_time)/60;
disp(['Total simulation time is ' ...
    num2str(total_time_in_minutes) ' minutes.']);
